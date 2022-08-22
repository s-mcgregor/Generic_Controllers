close all
clc 

% Define variables to be used throught problem
w = logspace(-3,4,1000);
dd=0.:.001:2*pi;
xx1=cos(dd)-1;yy1=sin(dd);% This creates a unit circle for the Nyquist plot

rtd = 180./pi;
grav = 32.2;

% Define System Parameters

Za_V = -1.3046;
Ma   = 47.711;
Zd_V = -.2142;
Md   = -104.83;
V    = 886.78;
Za   = V*Za_V;
Zd   = V*Zd_V;
Mq = 0;
w_act = 2.*pi*11;
z_act = 0.707;

% design a RSLQR command to track Az using state feeback controller
% Set up Plant Model 
% See slide 16 from the notes (Week 5)

Ap = [Za/V     Za      0                  Zd;
      Ma/Za     0     (Md-(Ma*Zd/Za))     0;
      0         0      0                  1;
      0         0      -w_act*w_act      -2*z_act*w_act];

Bp = [0.; 0.; 0.; w_act*w_act];

Cp = eye(4);

Dp = 0.*Cp*Bp;


% Form the wiggle system
% See slide 16 from the notes (Week 5)
Aw = [0     1         0      0                  0;
      0     Za/V     Za      0                  Zd;
      0     Ma/Za     0     (Md-(Ma*Zd/Za))     0;
      0     0         0      0                  1;
      0     0         0      -w_act*w_act      -2*z_act*w_act];

Bw = [0.; 0.; 0.; 0.; w_act*w_act ];
F = [-1; 0.; 0.; 0.; 0.];

% Setup range of penalties for the LQR
Q=0.*Aw; % sets up a matrix Q that is the size of Aw and is all 0s
R=1;

%qq is a vector which each element scales the LQR Q matrix
qq=logspace(-6,0,500); %these are the varying values of q11

%{

% The only penalty is the 1,1 element on the error state
% Design the gains

for ii = 1:length(qq)
    
    Q(1,1)=qq(ii);
    [Kx_lqr,~,~]=lqr(Aw,Bw,Q,R);

    % Form matrices for feedforward gain calculation
    A_p = Aw(2:5,2:5);
    B_p = Bw(2:5);
    C_p_reg = [Za 0 0 0];
    D_p_reg = Zd;
    Kfb = Kx_lqr(2:5);
    A_p_cl = A_p-B_p*Kfb;
    C_p_reg_cl = C_p_reg - D_p_reg*Kfb;
    Kp_ff = inv( -C_p_reg_cl*inv(A_p_cl)*B_p + D_p_reg );
    scale_Kff = 0.5; % Professor Wise recommended scaling down FF gain
    Kff = scale_Kff*Kp_ff;

    % populate the controller matrices
    Ac = 0.;
    Bc1 = [1. 0. 0. 0.];
    Bc2 = -1;
    Cc = -Kx_lqr(1);
    Dc1 = [-Kx_lqr(2:5)];
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
    dmax = max(abs(y(:,3)))*rtd*grav; % compute in per g commanded
    ddmax = max(abs(y(:,4)))*rtd*grav;
    metric=[qq(ii) rdmin srmin wc taur taus azmin azmax dmax ddmax Tmax Smax];
    data(ii,:) = metric;
    az_st(ii,:) = az';
    del_st(ii,:) = rtd*y(:,3);
    deldot_st(ii,:) = rtd*y(:,4);

end

% xlswrite('H5_Metrics.xls',data);
%}


% Set Q matrix based on selected value
Q(1,1)=qq(105);
[Kx_lqr,~,~]=lqr(Aw,Bw,Q,R);

% Form matrices for feedforward gain calculation
A_p = Aw(2:5,2:5);
B_p = Bw(2:5);
C_p_reg = [Za 0 0 0];
D_p_reg = Zd;
Kfb = Kx_lqr(2:5);
A_p_cl = A_p-B_p*Kfb;
C_p_reg_cl = C_p_reg - D_p_reg*Kfb;
Kp_ff = inv( -C_p_reg_cl*inv(A_p_cl)*B_p + D_p_reg );
scale_Kff = 0.5; % Professor Wise recommended scaling down FF gain
Kff = scale_Kff*Kp_ff;

% populate the controller matrices
Ac = 0.;
Bc1 = [1. 0. 0. 0.];
Bc2 = -1;
Cc = -Kx_lqr(1);
Dc1 = [-Kx_lqr(2:5)];
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

% SS model of loop gain at the plant input Lu
     A_Lu = [ Ap 0.*Bp*Cc;  Bc1*Cp Ac];
     B_Lu = [ Bp; Bc1*Dp];
     C_Lu = -[ Dc1*Cp Cc];%change sign for loop gain
     D_Lu = -[ Dc1*Dp];
     sys_Lu = ss(A_Lu,B_Lu,C_Lu,D_Lu);
%SS model of loop gain Ly at the plant output
     Aout = [ Ap Bp*Cc;  0.*Bc1*Cp Ac];
     Bout = [ Bp*Dc1; Bc1];
     Cout = -[ Cp Dp*Cc];%change sign for loop gain
     Dout = -[ Dp*Dc1];
     sys_Ly = ss(Aout,Bout,Cout,Dout);
% Analysis at plant input
     magdb = 20*log10(abs(squeeze(freqresp(sys_Lu,w))));
     wc = crosst(magdb,w); % LGCF, assumes Lu is a scalar
     sr = sigma(sys_Lu,w,3);
     sru_min = min(abs(sr));
     rd = sigma(sys_Lu,w,2);
     rdu_min = min(abs(rd));
     Lu = squeeze(freqresp(sys_Lu,w));
% Analysis at plant output
     T  = freqresp(sys_cl,w); % Complementary Sensitivity
     S = 1 - T; % Sensitivity
     T_Az = 20*log10(abs(squeeze(T(1,1,:))));
     S_Az = 20*log10(abs(squeeze(S(1,1,:))));
     Tmax = max(T_Az); % Inf Norm of T in dB
     Smax = max(S_Az); % Inf Norm of S in dB
     Ws_magdb = 20*log10(abs(squeeze(1/freqresp(Ws.sys,w)))) ;
     Wt_magdb = 20*log10(abs(squeeze(1/freqresp(Wt.sys,w)))) ;

     
disp('Homework 5, 5-State RSLQR')
    %Compute singluar value margins
     neg_gm =  min([ (1/(1+rdu_min)) (1-sru_min)]); % in dB
     pos_gm =  max([ (1/(1-rdu_min)) (1+sru_min)]); % in dB
     neg_gmdB = 20*log10( neg_gm ); % in dB
     pos_gmdB = 20*log10( pos_gm ); % in dB
     pm = 180*(max([2*asin(rdu_min/2) 2*asin(sru_min/2)]))/pi;% in deg

     disp('Singular value margins')
     disp(['Min Singular value I+Lu =    ' num2str(rdu_min)])
     disp(['Min Singular value I+invLu = ' num2str(sru_min)])
     disp(['Singular value gain margins = [' ...
           num2str(neg_gmdB) ' dB,' num2str(pos_gmdB) ' dB ]' ])
     disp(['Singular value phase margins = [ +/-' ...
           num2str(pm)  ' deg ]' ])     
disp(' ')      
     
% Time Domain Analysis
t = [0:.01:2];   %Time Vector

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
dmax = max(abs(y(:,3)))*rtd*grav; % compute in per g commanded
ddmax = max(abs(y(:,4)))*rtd*grav;

FIVES.y      = y;
FIVES.t      = t;
FIVES.az     = y(:,1); %  acceleration (fps2)
FIVES.taur   = taur;
FIVES.taus   = taus;
FIVES.Lu     = Lu;
FIVES.rd     = rd;
FIVES.rdu_min   = rdu_min;
FIVES.sr     = sr;
FIVES.sru_min   = sru_min;
FIVES.S      = S;
FIVES.T      = T;
FIVES.lgcf   = wc;
FIVES.dele_Az  = dele_Az;
FIVES.dele_q   = dele_q;
FIVES.deledot_Az  = deledot_Az;
FIVES.deledot_q   = deledot_q;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                         SPC                    %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Full State Feedback')
[V,D] = eig(Aw - Bw*Kx_lqr);
disp('Closed Loop Eigenvectors')
V
disp('Closed Loop Eigenvalues')
diag(D)

% Design SPC Controller

Eig_Keep = D(3:5,3:5); % keep dominant eigenvalues
X_r = V(:,3:5); % keep associated eigenvectors
C = [1 0 0 0 0;
     0 1 0 0 0;
     0 0 1 0 0];

K_bar = real(Kx_lqr*X_r*inv(C*X_r)); % Ky

% populate the controller matrices
% two extra zeros added on Bc1 and Dc1 to connect to 4 state plant
Ac = 0.;
Bc1 = [1. 0. 0. 0.];
Bc2 = -1.;
Cc = -K_bar(1);
Dc1 = [-K_bar(2:3) 0. 0.];
Dc2 = 0.;

Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
Acl = [ (Ap+Bp*Z*Dc1*Cp) (Bp*Z*Cc)
    (Bc1*(Cp+Dp*Z*Dc1*Cp)) (Ac+Bc1*Dp*Z*Cc)];
Bcl = [ Bp*Z*Dc2;
    (Bc2+Bc1*Dp*Z*Dc2)];
Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
Dcl =(Dp*Z*Dc2);
sys_cl = ss(Acl,Bcl,Ccl,Dcl);

disp('Output Feedback')
[V,D] = eig(Aw - Bw*K_bar*C);
disp('Closed Loop Eigenvectors')
V
disp('Closed Loop Eigenvalues')
diag(D)


% SS model of loop gain at the plant input Lu
     A_Lu = [ Ap 0.*Bp*Cc;  Bc1*Cp Ac];
     B_Lu = [ Bp; Bc1*Dp];
     C_Lu = -[ Dc1*Cp Cc];%change sign for loop gain
     D_Lu = -[ Dc1*Dp];
     sys_Lu = ss(A_Lu,B_Lu,C_Lu,D_Lu);
%SS model of loop gain Ly at the plant output
     Aout = [ Ap Bp*Cc;  0.*Bc1*Cp Ac];
     Bout = [ Bp*Dc1; Bc1];
     Cout = -[ Cp Dp*Cc];%change sign for loop gain
     Dout = -[ Dp*Dc1];
     sys_Ly = ss(Aout,Bout,Cout,Dout);
% Analysis at plant input
     magdb = 20*log10(abs(squeeze(freqresp(sys_Lu,w))));
     wc = crosst(magdb,w); % LGCF, assumes Lu is a scalar
     sr = sigma(sys_Lu,w,3);
     sru_min = min(abs(sr));
     rd = sigma(sys_Lu,w,2);
     rdu_min = min(abs(rd));
     Lu = squeeze(freqresp(sys_Lu,w));
% Analysis at plant output
     T  = freqresp(sys_cl,w); % Complementary Sensitivity
     S = 1 - T; % Sensitivity
     T_Az = 20*log10(abs(squeeze(T(1,1,:))));
     S_Az = 20*log10(abs(squeeze(S(1,1,:))));
     Tmax = max(T_Az); % Inf Norm of T in dB
     Smax = max(S_Az); % Inf Norm of S in dB
     Ws_magdb = 20*log10(abs(squeeze(1/freqresp(Ws.sys,w)))) ;
     Wt_magdb = 20*log10(abs(squeeze(1/freqresp(Wt.sys,w)))) ;     
     
     
disp('Homework 5, SPC')
    %Compute singluar value margins
     neg_gm =  min([ (1/(1+rdu_min)) (1-sru_min)]); % in dB
     pos_gm =  max([ (1/(1-rdu_min)) (1+sru_min)]); % in dB
     neg_gmdB = 20*log10( neg_gm ); % in dB
     pos_gmdB = 20*log10( pos_gm ); % in dB
     pm = 180*(max([2*asin(rdu_min/2) 2*asin(sru_min/2)]))/pi;% in deg

     disp('Singular value margins')
     disp(['Min Singular value I+Lu =    ' num2str(rdu_min)])
     disp(['Min Singular value I+invLu = ' num2str(sru_min)])
     disp(['Singular value gain margins = [' ...
           num2str(neg_gmdB) ' dB,' num2str(pos_gmdB) ' dB ]' ])
     disp(['Singular value phase margins = [ +/-' ...
           num2str(pm)  ' deg ]' ])         
disp(' ')        
     
% Time Domain Analysis
t = [0:.01:2];   %Time Vector

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
dmax = max(abs(y(:,3)))*rtd*grav; % compute in per g commanded
ddmax = max(abs(y(:,4)))*rtd*grav;

% Compute the noise-to-control TF
 Bv = [       Bp*Z*Dc1;
     (Bc1+Bc1*Dp*Z*Dc1)];
 Cv  = [ Z*Dc1*Cp Z*Cc];
 Cvv = [ Cv ; Cv*Acl ];
 Dv = Z*Dc1;
 Dvv = [ Dv; Cv*Bv];
 sys_noise = ss(Acl,Bv,Cvv,Dvv);
 v_2_u  = freqresp(sys_noise,w); % Noise to control freq response
 dele_Az    = 20*log10(abs(squeeze(v_2_u(1,1,:))));
 dele_q     = 20*log10(abs(squeeze(v_2_u(1,3,:))));
 deledot_Az = 20*log10(abs(squeeze(v_2_u(2,1,:))));
 deledot_q  = 20*log10(abs(squeeze(v_2_u(2,3,:))));


SPC.y      = y;
SPC.t      = t;
SPC.az     = y(:,1); %  acceleration (fps2)
SPC.taur   = taur;
SPC.taus   = taus;
SPC.Lu     = Lu;
SPC.rd     = rd;
SPC.rdu_min   = rdu_min;
SPC.sr     = sr;
SPC.sru_min   = sru_min;
SPC.S      = S;
SPC.T      = T;
SPC.lgcf   = wc;
SPC.dele_Az  = dele_Az;
SPC.dele_q   = dele_q;
SPC.deledot_Az  = deledot_Az;
SPC.deledot_q   = deledot_q;

%{
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

%}

% Az comparison
figure('Name','Acceleration Time History')
plot(Hmk3_RSLQR.t, Hmk3_RSLQR.y(:,1))
hold on
plot(HINF.t,HINF.az);
plot(t,FIVES.az);
plot(t,SPC.az);
xlabel('Time (sec)');
ylabel('Az (fps2)');
legend(['3 State RSLQR 63% Tr = ' num2str(Hmk3_RSLQR.metric(4)) ' 95% Ts = ' num2str(Hmk3_RSLQR.metric(5))],...
       ['Hinf 63% Tr = ' num2str(HINF.taur) ' 95% Ts = ' num2str(HINF.taus)],...
       ['5 State RSLQR 63% Tr = ' num2str(FIVES.taur) ' 95% Ts = ' num2str(FIVES.taus)],...
       ['SPC 63% Tr = ' num2str(SPC.taur) ' 95% Ts = ' num2str(SPC.taus)],...
       'Location','Best');
hold off


% Plot Nyquist of Lu
% Need to define a unit circle for nyquist plot
n = 36; %// Define number of points on circle
theta = linspace(0, 2*pi, n);
x = cos(theta);
y = sin(theta);
figure('Name','Nyquist Plot at Plant Input'),
plot(x-1, y, '-.r'); % Unit circle
hold on
plot(real(squeeze(Hmk3_RSLQR.Lu)),imag(squeeze(Hmk3_RSLQR.Lu)));
plot(real(squeeze(HINF.Lu)),imag(squeeze(HINF.Lu)));
plot(real(squeeze(FIVES.Lu)),imag(squeeze(FIVES.Lu)));
plot(real(squeeze(SPC.Lu)),imag(squeeze(SPC.Lu)));
hold off
axis([-2 2 -2 2]);
legend('Unit Circle','3 State RSLQR','Hinf-SF','5 State RSLQR','SPC','Location','Best');
xlabel('Re(L)')
ylabel('Im(L)')
title('Nyquist Plot at Plant Input')



figure('Name','Bode Magnitude at Plant Input'),
semilogx(w,20*log10(abs(Hmk3_RSLQR.Lu)),...
    w,20*log10(abs(HINF.Lu)),...
    w,20*log10(abs(FIVES.Lu)),...
    w,20*log10(abs(SPC.Lu)));
legend([ '3 State RSLQR LGCF = '  num2str(Hmk3_RSLQR.metric(3)) ],...
       [ 'Hinf-SF LGCF = '     num2str(HINF.metric(4)) ],...
       [ '5 State RSLQR LGCF = '  num2str(FIVES.lgcf)],...
       [ 'SPC LGCF = '  num2str(SPC.lgcf) ],...
       'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag')
title('Bode at Input')


figure('Name','Return Difference at Plant Input'),
semilogx(w,20*log10(abs(Hmk3_RSLQR.rd)),...
    w,20*log10(abs(HINF.rd)),...
    w,20*log10(abs(FIVES.rd)),...
    w,20*log10(abs(SPC.rd)))
legend(['3 State RSLQR min(I+Lu) = ' num2str(Hmk3_RSLQR.metric(1)) ],...
       [' Hinf min(I+Lu) = ' num2str(HINF.rdu_min) ],...
       ['5 State RSLQR min(I+Lu) = ' num2str(FIVES.rdu_min) ],...
       ['SPC min(I+Lu) = ' num2str(SPC.rdu_min)],...
       'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Return Difference at Plant Input')


figure('Name','Stability Robustness at Plant Input'),
semilogx(w,20*log10(abs(Hmk3_RSLQR.sr)),...
    w,20*log10(abs(HINF.sr)),...
    w,20*log10(abs(FIVES.sr)),...
    w,20*log10(abs(SPC.sr)))
legend(['3 State RSLQR min(I+invLu) = ' num2str(Hmk3_RSLQR.metric(2)) ],...
       [' Hinf min(I+invLu) = ' num2str(HINF.sru_min) ],...
       ['5 State RSLQR min(I+invLu) = ' num2str(FIVES.sru_min) ],...
       ['SPC min(I+invLu) = ' num2str(SPC.sru_min)],...
       'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Stability Robustness at Plant Input')

% Complementary Sensitivity plot
figure('Name','Comp Sensitivity T=Y/R')
T_az_max = max(abs(squeeze(HINF.T(1,1,:))));
semilogx(w,20*log10(abs(squeeze(HINF.T(1,1,:)))),'LineWidth',2);grid
tex1 = ['T-Az inf norm = ' num2str(T_az_max)];
hold on
T_az_max_RSLQR = max(abs(squeeze(Hmk3_RSLQR.T(1,1,:))));
semilogx(w,20*log10(abs(squeeze(Hmk3_RSLQR.T(1,1,:)))),'LineWidth',2);grid
tex2 = ['T-Az norm 3 State RSLQR= ' num2str(T_az_max_RSLQR)];

T_az_max_Fives = max(abs(squeeze(FIVES.T(1,1,:))));
semilogx(w,20*log10(abs(squeeze(FIVES.T(1,1,:)))),'LineWidth',2);grid
tex3 = ['T-Az norm 5 State RSLQR= ' num2str(T_az_max_Fives)];

T_az_max_SPC = max(abs(squeeze(SPC.T(1,1,:))));
semilogx(w,20*log10(abs(squeeze(SPC.T(1,1,:)))),'LineWidth',2);grid
tex4 = ['T-Az norm SPC= ' num2str(T_az_max_SPC)];

legend(tex1, tex2, tex3, tex4, 'Location','Best');
xlabel('Frequency (rps)')
ylabel('Co-Sens Mag (dB)')
title('Comp Sensitivity T=Y/R')
hold off

% Sensitivity plot
figure('Name',' Sensitivity S=E/R')
S_az_max = max(abs(squeeze(HINF.S(1,1,:))));
semilogx(w,20*log10(abs(squeeze(HINF.S(1,1,:)))),'LineWidth',2);grid
tex1 = ['S-Az INF norm = ' num2str(S_az_max)];
hold on
S_az_max_RSLQR = max(abs(squeeze(Hmk3_RSLQR.S(1,1,:))));
semilogx(w,20*log10(abs(squeeze(Hmk3_RSLQR.S(1,1,:)))),'LineWidth',2);grid
tex2 = ['S-Az norm 3 State RSLQR= ' num2str(S_az_max_RSLQR)];

S_az_max_FIVES = max(abs(squeeze(FIVES.S(1,1,:))));
semilogx(w,20*log10(abs(squeeze(FIVES.S(1,1,:)))),'LineWidth',2);grid
tex3 = ['S-Az norm 5 State RSLQR= ' num2str(S_az_max_FIVES)];

S_az_max_SPC = max(abs(squeeze(SPC.S(1,1,:))));
semilogx(w,20*log10(abs(squeeze(SPC.S(1,1,:)))),'LineWidth',2);grid
tex4 = ['S-Az norm SPC= ' num2str(S_az_max_SPC)];

legend(tex1, tex2, tex3, tex4, 'Location','Best');
xlabel('Frequency (rps)')
ylabel('Co-Sens Mag (dB)')
title('Sensitivity S=E/R')
hold off

% Plot Noise-To-Control
figure('Name','Noise-2-Control');
semilogx(w,HINF.dele_Az,'k',w,HINF.dele_q,'k--',...
         w,Hmk3_RSLQR.dele_Az,'b',w,Hmk3_RSLQR.dele_q,'b--',...
         w,FIVES.dele_Az,'g',w,FIVES.dele_q,'g--',...
         w,SPC.dele_Az,'m',w,SPC.dele_q,'m--',...
         'LineWidth',2);
legend('Az Hinf','q Hinf',...
       'Az RSlqr 3 State','q RSlqr 3 State',...
       'Az RSlqr 5 State','q RSlqr 5 State',...
       'Az SPC','q SPC',...
       'Location','Best');
title(['Noise to Control ']);
xlabel('Freq (rps)');ylabel('Mag (dB)');
grid;


% Plot Noise-To-Control Rate
figure('Name','Noise-2-Control Rate');
semilogx(w,HINF.deledot_Az,'k',w,HINF.deledot_q,'k--',...
         w,Hmk3_RSLQR.deledot_Az,'b',w,Hmk3_RSLQR.deledot_q,'b--',...
         w,FIVES.deledot_Az,'g',w,FIVES.deledot_q,'g--',...
         w,SPC.deledot_Az,'m',w,SPC.deledot_q,'m--',...
         'LineWidth',2);
legend('Az Hinf','q Hinf',...
       'Az RSlqr 3 State','q RSlqr 3 State',...
       'Az RSlqr 5 State','q RSlqr 5 State',...
       'Az SPC','q SPC',...
       'Location','Best');
title(['Noise to Control Rate']);
xlabel('Freq (rps)');ylabel('Mag (dB)');
grid;