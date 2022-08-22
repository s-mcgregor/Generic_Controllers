clear all
close all

% Define variables to be used throught problem
w = logspace(-2,3,500);
dd=0.:.001:2*pi;
xx1=cos(dd)-1;yy1=sin(dd);% This creates a unit circle for teh Nyquist plot

%t = 0.:0.01:5.;
rtd = 180./pi;

% ME598 Robust Control - Homework 3, Problem 1

Za_V = -1.3046;
Ma   = 47.711;
Zd_V = -.2142;
Md   = -104.83;
V    = 886.78;
Za   = V*Za_V;
Zd   = V*Zd_V;
Mq = -1.03341;
w_act = 2.*pi*11;
z_act = 0.707;      

Ap = [Za_V  1.      Zd_V          0.;

      Ma    0.       Md           0.;

      0.    0.       0.           1.;

      0.    0. -w_act*w_act -2*z_act*w_act];

Bp = [0.; 0.; 0.; w_act*w_act ]; 

Cp = [  Za   0.  Zd  0.;

        eye(4)];

Dp = 0.*Cp*Bp;

[~,nu] = size(Bp);
[ny,nx] = size(Cp);

% Compute open loop eigenvalues to determine if plant is unstable
disp('   ');
disp('Open Loop Eigenvalues')
damp(Ap)

disp('   ')
disp('Plant TFr')
sys_p = ss(Ap,Bp,Cp,Dp);
zpk(sys_p)
disp('   ')

for jj = 1:nu,
for ii = 1:ny,
    disp(['Finitie Zeros for Input = ' num2str(jj) ' Output = ' num2str(ii)])
    tzero(Ap,Bp(:,jj),Cp(ii,:),Dp(ii,jj))
end
end

% design a RSLQR command to track Az using state feeback controller
% form the wiggle system:
% state vector: [int(e) (fps), AOA (rad), pitch rate q (rps)]'
Aw = [0. Za 0.;
0. Za_V 1.;
0. Ma 0.];
Bw = [ Zd; Zd_V; Md];

% Setup range of penalties for the LQR
Q=0.*Aw; %sets up a matrix Q that is the size of Aw and is all 0s
R=1;

%qq is a vector which each element scales the LQR Q matrix
qq=logspace(-6,0,500); %these are the varying values of q11

% The only penalty is the 1,1 element on the error state
% Design the gains
for ii = 1:length(qq)
    
    Q(1,1)=qq(ii);
    [Kx_lqr,~,~]=lqr(Aw,Bw,Q,R);

    % Form matrices for feedforward gain calculation
    A_p = Aw(2:3,2:3);
    B_p = Bw(2:3);
    C_p_reg = [Za 0];
    D_p_reg = Zd;
    Kfb = Kx_lqr(2:3);
    A_p_cl = A_p-B_p*Kfb;
    C_p_reg_cl = C_p_reg - D_p_reg*Kfb;
    Kp_ff = inv( -C_p_reg_cl*inv(A_p_cl)*B_p + D_p_reg );
    scale_Kff = 0.5; % Professor Wise recommended scaling down FF gain
    Kff = scale_Kff*Kp_ff;

    % populate the controller matrices
    Ac = 0.;
    Bc1 = [1. 0. 0. 0. 0.];
    Bc2 = -1;
    Cc = -Kx_lqr(1);
    Dc1 = [0. -Kx_lqr(2:3) 0. 0.];
    Dc2 = Kff;

    % Form the closed loop system
    Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
    Acl = [ (Ap+Bp*Z*Dc1*Cp) (Bp*Z*Cc);
    (Bc1*(Cp+Dp*Z*Dc1*Cp)) (Ac+Bc1*Dp*Z*Cc)];
    Bcl = [ Bp*Z*Dc2;
    (Bc2+Bc1*Dp*Z*Dc2)];
    Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
    Dcl =(Dp*Z*Dc2);
    sys_cl = ss(Acl,Bcl,Ccl,Dcl);

    % Frequency Domain Analysis
    % SS model of loop gain at the plant input Lu
    A_Lu = [ Ap 0.*Bp*Cc; Bc1*Cp Ac];
    B_Lu = [ Bp; Bc1*Dp];
    C_Lu = -[ Dc1*Cp Cc];%change sign for loop gain
    D_Lu = -[ Dc1*Dp];
    sys_Lu = ss(A_Lu,B_Lu,C_Lu,D_Lu);

    t = [0:.01:2];   %Time Vector

    magdb = 20*log10(abs(squeeze(freqresp(sys_Lu,w))));
    wc = crosst(magdb,w); % LGCF, assumes Lu is a scalar
    sr = sigma(sys_Lu,w,3); % Stability Robustness
    srmin = min(abs(sr));
    rd = sigma(sys_Lu,w,2); % Return Difference
    rdmin = min(abs(rd));
    T = freqresp(sys_cl,w); % Complementary Sensitivity
    S = 1 - T; % Sensitivity
    T_st(ii,:) = 20*log10(abs(squeeze(T(1,1,:))));
    S_st(ii,:) = 20*log10(abs(squeeze(S(1,1,:))));
    Tmax = max(T_st(ii,:)); % Inf Norm of T in dB
    Smax = max(S_st(ii,:)); % Inf Norm of S in dB

    % Time Domain Analysis
    y = step(sys_cl,t);
    az = y(:,1); % acceleration (fps2)
    aze = abs(ones(size(az))-az); % error for az
    taur = 0.; taus= 0.; % rise time and settling time
    fv = aze(numel(aze)); % final value of the error
    e_n = aze - fv*ones(size(aze)) - 0.36*ones(size(aze));
    e_n1 = abs(e_n) + e_n;
    taur = crosst(e_n1,t); % rise time
    e_n = aze - fv*ones(size(aze)) - 0.05*ones(size(aze));
    e_n1 = abs(e_n) + e_n;
    taus = crosst(e_n1,t); % settling time
    azmin = abs(min(az))*100; % undershoot
    azmax = (abs(max(az))-1)*100; % overshoot
    grav = 32.2; %fps2
    dmax = max(abs(y(:,4)))*rtd*grav; % compute in per g commanded
    ddmax = max(abs(y(:,5)))*rtd*grav;
    metric=[qq(ii) rdmin srmin wc taur taus azmin azmax dmax ddmax Tmax Smax];
    data(ii,:) = metric;
    az_st(ii,:) = az';
    del_st(ii,:) = rtd*y(:,4);
    deldot_st(ii,:) = rtd*y(:,5);

end

Hmk3_RSLQR.t = t;
Hwk3_RSLQR.y = y;

%for ii = 1:length(qq)
for ii = 1:165
    plot(t,az_st(ii,:))
    if ii == 1
       hold on
    end
end
xlabel('Time (s)')
ylabel('Az (fps2)')
hold off

%for ii = 1:length(qq)
figure
for ii = 1:165
    plot(t,del_st(ii,:))
    if ii == 1
       hold on
    end
end
xlabel('Time (s)')
ylabel('Elevon (deg)')
hold off

%for ii = 1:length(qq)
figure
for ii = 1:165
    plot(t,deldot_st(ii,:))
    if ii == 1
       hold on
    end
end
xlabel('Time (s)')
ylabel('Elevon Rate (dps)')
hold off

% Only track at selected qq [qq(171) is original selection] 165 best
% selection... 234 best performance
Q=0.*Aw; %sets up a matrix Q that is the size of Aw and is all 0s
Q(1,1)=qq(165);
[Kx_lqr,~,~]=lqr(Aw,Bw,Q,R);
Kx_lqr
% Form matrices for feedforward gain calculation
A_p = Aw(2:3,2:3);
B_p = Bw(2:3);
C_p_reg = [Za 0];
D_p_reg = Zd;
Kfb = Kx_lqr(2:3);
A_p_cl = A_p-B_p*Kfb;
C_p_reg_cl = C_p_reg - D_p_reg*Kfb;
Kp_ff = inv( -C_p_reg_cl*inv(A_p_cl)*B_p + D_p_reg );
scale_Kff = .5; % Professor Wise recommended scaling down FF gain
Kff = scale_Kff*Kp_ff;

% populate the controller matrices
Ac = 0.
Bc1 = [1. 0. 0. 0. 0.]
Bc2 = -1
Cc = -Kx_lqr(1)
Dc1 = [0. -Kx_lqr(2:3) 0. 0.]
Dc2 = Kff

% Form the closed loop system
Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
Acl = [ (Ap+Bp*Z*Dc1*Cp) (Bp*Z*Cc);
    (Bc1*(Cp+Dp*Z*Dc1*Cp)) (Ac+Bc1*Dp*Z*Cc)];
Bcl = [ Bp*Z*Dc2;
    (Bc2+Bc1*Dp*Z*Dc2)];
Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
Dcl =(Dp*Z*Dc2);
    sys_cl = ss(Acl,Bcl,Ccl,Dcl);

% Frequency Domain Analysis
% SS model of loop gain at the plant input Lu
A_Lu = [ Ap 0.*Bp*Cc; Bc1*Cp Ac];
B_Lu = [ Bp; Bc1*Dp];
C_Lu = -[ Dc1*Cp Cc];%change sign for loop gain
D_Lu = -[ Dc1*Dp];

L_input = ss(A_Lu,B_Lu,C_Lu,D_Lu);
S_input = inv(eye(size(L_input))+L_input); % Sensitivity
T_input = eye(size(L_input)) - S_input; % Complimentary Sensitivity
Neg_T_input =-1*T_input;

figure
sigma(inv(Neg_T_input))
hold on

wn = 69.11; % 11 Hz; 69.11 rps
damp = .707;
num = [-1 -2*damp*wn 0];
den = [1 2*damp*wn wn*wn];
D01 = tf(num,den);
sigma(D01)
hold off
legend('1/sigma(M)','sigma(Delta)')

% Break loop at plant input (Lu); Plot Bode Plot, id gain and phase margins, LGXF, PXF
figure
margin(L_input)
[GM, PM_deg, wc_GM, wc_Pm] = margin(L_input);

% Plot Nyquist of Lu, id gain and phase margin, loop gain and phase xover
%freq
% Need to define a unit circle for nyquist plot
n = 36; %// Define number of points on circle
theta = linspace(0, 2*pi, n);
x = cos(theta);
y = sin(theta);

figure
nyquist(L_input)
grid on
axis([-3 3 -3 3])
hold on;
plot(x-1, y, '-.r'); % Unit circle
hold off;

%Sigma(I+Lu)
figure
sigma(eye(size(L_input))+L_input);
title('Sigma (I + Lu)')
%legend([' min(I+Lu) = ' num2str(RDu_min)],'Location','Best');

%Sigma(I+Inv(Lu))
figure
sigma(eye(size(L_input))+inv(L_input));
title('Sigma (I + inv(Lu))')
%legend([' min(I+invLu) = ' num2str(SRu_min)],'Location','Best');

% Sensitivity and Complimentary Sensitivity plots
HT = freqresp(T_input, w);
HT_squeeze = squeeze(HT); %squeeze function reduces dimensions
HS = freqresp(S_input, w);
HS_squeeze = squeeze(HS); %squeeze function reduces dimensions

figure
plot(w,abs(HT_squeeze))
hold on
plot(w,abs(HS_squeeze))
xlim([0 250])
xlabel('Frequency (rad/s)')
ylabel('Response')
title('Sensitivity (S) and Complimentary Sensitivity (T)')
legend('T','S')
hold off

Cw = eye(size(Aw)); %output all states
Dw = 0.*Cw*Bw;
figure
rlocus(ss(Aw,Bw,Cw,Dw))