clear all
close all

% Define variables to be used throught problem
w = logspace(-3,4,1000);
dd=0.:.001:2*pi;
xx1=cos(dd)-1;yy1=sin(dd);% This creates a unit circle for the Nyquist plot

%t = 0.:0.01:5.;
rtd = 180./pi;
grav = 32.2;

% ME598 Robust Control - Homework 3, Problem 1

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

Aw = [0 1 0 0 0;
     0 Za_V  1.      Zd_V          0.;

     0 Ma    0.       Md           0.;

     0 0.    0.       0.           1.;

     0 0.    0. -w_act*w_act -2*z_act*w_act];

Bw = [0.; 0.; 0.; 0.; w_act*w_act ]; 

Q=0.*Aw; %sets up a matrix Q that is the size of Aw and is all 0s
Q(1,1)=.00009374;
R=1;
Kx_lqr = [ 0.0097 -2.6491 -0.2393];

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

%**************************************************************************
% Compute Frequency domain results for LQR and LQG using looping method
%**************************************************************************
% Close the loop
% Plant form  xdot = Apx + Bpu;
%                      y = Cpx +Dpu
% Controller xcdot = Acxc + Bc1y + Bc2r
%                      u = Ccxc + Dc1y + Dc2r

Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
Acl = [     (Ap+Bp*Z*Dc1*Cp)              (Bp*Z*Cc);
    (Bc1*(Cp+Dp*Z*Dc1*Cp))  (Ac+Bc1*Dp*Z*Cc)];
Bcl = [       Bp*Z*Dc2;
    (Bc2+Bc1*Dp*Z*Dc2)];
Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
Dcl =(Dp*Z*Dc2);
Ccl_Az = Ccl(1,:);
Dcl_Az =Dcl(1,:);
sys_cl = ss(Acl,Bcl,Ccl,Dcl);

% % Time Domain Analysis
t = [0:.01:2];   %Time Vector
[y_cl,x_cl] = step(Acl,Bcl,Ccl_Az,Dcl_Az,1,t);
sys_cl      = ss(Acl,Bcl,Ccl_Az,Dcl_Az);

% Preallocate for speed
Lu     = 0.*w;
RDu    = 0.*w;
SRu    = 0.*w;
RDy2_min = 0.*w;
Sy_1   = 0.*w;
Ty_1   = 0.*w;
Sy_2   = 0.*w;
Ty_2   = 0.*w;

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


% Analysis at plant inut
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

y = step(sys_cl,t);
az = y(:,1); %  acceleration (fps2)
aze = abs(ones(size(az))-az);  % error for az
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
%dmax = max(abs(y(:,3)))*rtd*grav; % compute in per g commanded
%ddmax = max(abs(y(:,4)))*rtd*grav;
metric=[rdu_min sru_min wc taur taus azmin azmax Tmax Smax];

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

Hmk3_RSLQR.y      = y;
Hmk3_RSLQR.t      = t;
Hmk3_RSLQR.Lu     = Lu;
Hmk3_RSLQR.rd     = rd;
Hmk3_RSLQR.sr     = sr;
Hmk3_RSLQR.S      = S;
Hmk3_RSLQR.T      = T;
Hmk3_RSLQR.dele_Az  = dele_Az;
Hmk3_RSLQR.dele_q   = dele_q;
Hmk3_RSLQR.deledot_Az  = deledot_Az;
Hmk3_RSLQR.deledot_q   = deledot_q;
Hmk3_RSLQR.metric = metric;     