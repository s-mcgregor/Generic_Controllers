close all
clear all
clc 

% Define variables to be used throught problem
w = logspace(-2,3,500);
dd=0.:.001:2*pi;
xx1=cos(dd)-1;yy1=sin(dd);% This creates a unit circle for the Nyquist plot

%t = 0.:0.01:5.;
rtd = 180./pi;

% Set Parameters

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

% design a RSLQR command to track Az using state feeback controller
% form the wiggle system:
% state vector: [int(e) (fps), AOA (rad), pitch rate q (rps)]'
Aw = [0. Za 0.;
0. Za_V 1.;
0. Ma 0.];
Bw = [ Zd; Zd_V; Md];

% Setup range of penalties for the LQR
Q_LQR=0.*Aw; %sets up a matrix Q that is the size of Aw and is all 0s
Q_LQR(1,1)=.00009374;
R_LQR=1;
% Solve for RSLQR state feedback gains
[Kc,Pw]=lqr(Aw,Bw,Q_LQR,R_LQR);

format short G
disp('LQR Gains');
Kc

% Form matrices for feedforward gain calculation
A_p = Aw(2:3,2:3);
B_p = Bw(2:3);
C_p_reg = [Za 0];
D_p_reg = Zd;
Kfb = Kc(2:3);
A_p_cl = A_p-B_p*Kfb;
C_p_reg_cl = C_p_reg - D_p_reg*Kfb;
Kp_ff = inv( -C_p_reg_cl*inv(A_p_cl)*B_p + D_p_reg );
scale_Kff = .5; % Professor Wise recommended scaling down FF gain
Kff = scale_Kff*Kp_ff;

% populate the controller matrices
Ac = 0.;
Bc1 = [1. 0. 0. 0. 0.];
Bc2 = -1;
Cc = -Kc(1);
Dc1 = [0. -Kc(2:3) 0. 0.];
Dc2 = Kff;

% Connect the controller to the plant model

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
     
disp('RSLQR margins')
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
q = y(:,2); % pitch rate (dps)
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

RSLQR.y      = y;
RSLQR.t      = t;
RSLQR.az     = y(:,1); %  acceleration (fps2)
RSLQR.q      = y(:,2); %  pitch rate (dps)
RSLQR.dele   = y(:,3); %  dele (deg)
RSLQR.dd     = y(:,4); %  dele_dot (dps)
RSLQR.taur   = taur;
RSLQR.taus   = taus;
RSLQR.Lu     = Lu;
RSLQR.rd     = rd;
RSLQR.rdu_min   = rdu_min;
RSLQR.sr     = sr;
RSLQR.sru_min   = sru_min;
RSLQR.S      = S;
RSLQR.T      = T;
RSLQR.lgcf   = wc;
RSLQR.dele_Az  = dele_Az;
RSLQR.dele_q   = dele_q;
RSLQR.deledot_Az  = deledot_Az;
RSLQR.deledot_q   = deledot_q;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Next is the LQG design. The nominal process and measurement noise
% covariance matrices are given in the homework

C = [1 0 0; 0 0 1];
G = eye(3);

Q0 = diag([0.001 0.0014 0.005]);
R0 = 1000*diag([0.025^2 0.001^2]);
rho = [10^10 10^4 10^2];


Qf_10 = Q0 + (1/rho(1))*Bw*Bw';
Qf_4 = Q0 + (1/rho(2))*Bw*Bw';
Qf_2 = Q0 + (1/rho(3))*Bw*Bw';

% Q for rho = 10^10
[L, Pf] = lqe(Aw,G,C,Qf_10,R0);
disp('Rho = 10^10');

% L2 = (lqr(Aw',C',Qf_10,R0))' <-- This is another way of writing it with
% the lqr command

Ac = [ 0 0 0 0;
    L(:,1) (Aw - Bw*Kc - L*C)];

Bc1 = [1 0 0 0 0;
    zeros(3,2) L(:,2) zeros(3,2)];

Bc2 = [-1; -1; 0 ; 0];

Cc = [ 0 -Kc];

Dc1 = zeros(1,5);

Dc2 = 0;

% To form a closed loop system, the controller is connected to the plant

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
     
disp('LQE; rho = 10^10 margins');
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

LQE10.y      = y;
LQE10.t      = t;
LQE10.az     = y(:,1); %  acceleration (fps2)
LQE10.q      = y(:,2); %  pitch rate (dps)
LQE10.dele   = y(:,3); %  dele (deg)
LQE10.dd     = y(:,4); %  dele_dot (dps)
LQE10.taur   = taur;
LQE10.taus   = taus;
LQE10.Lu     = Lu;
LQE10.rd     = rd;
LQE10.rdu_min   = rdu_min;
LQE10.sr     = sr;
LQE10.sru_min   = sru_min;
LQE10.S      = S;
LQE10.T      = T;
LQE10.lgcf   = wc;
LQE10.dele_Az  = dele_Az;
LQE10.dele_q   = dele_q;
LQE10.deledot_Az  = deledot_Az;
LQE10.deledot_q   = deledot_q;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Q for rho = 10^4
[L, Pf] = lqe(Aw,G,C,Qf_4,R0);
disp('Rho = 10^4');
disp('Kalman Filter Gains');
L;

% L2 = (lqr(Aw',C',Qf_10,R0))' <-- This is another way of writing it with
% the lqr command

Ac = [ 0 0 0 0;
    L(:,1) (Aw - Bw*Kc - L*C)];

Bc1 = [1 0 0 0 0;
    zeros(3,2) L(:,2) zeros(3,2)];

Bc2 = [-1; -1; 0 ; 0];

Cc = [ 0 -Kc];

Dc1 = zeros(1,5);

Dc2 = 0;

% To form a closed loop system, the controller is connected to the plant

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
     
disp('LQE; rho = 10^4 margins')
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

LQE4.y      = y;
LQE4.t      = t;
LQE4.az     = y(:,1); %  acceleration (fps2)
LQE4.q      = y(:,2); %  pitch rate (dps)
LQE4.dele   = y(:,3); %  dele (deg)
LQE4.dd     = y(:,4); %  dele_dot (dps)
LQE4.taur   = taur;
LQE4.taus   = taus;
LQE4.Lu     = Lu;
LQE4.rd     = rd;
LQE4.rdu_min   = rdu_min;
LQE4.sr     = sr;
LQE4.sru_min   = sru_min;
LQE4.S      = S;
LQE4.T      = T;
LQE4.lgcf   = wc;
LQE4.dele_Az  = dele_Az;
LQE4.dele_q   = dele_q;
LQE4.deledot_Az  = deledot_Az;
LQE4.deledot_q   = deledot_q;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Q for rho = 10^2
[L, Pf] = lqe(Aw,G,C,Qf_2,R0);
disp('Rho = 10^2');
disp('Kalman Filter Gains');
L
disp('Kalman Filter SS Covariance');
Pf
Aw
C
Qf_2
R0

% L2 = (lqr(Aw',C',Qf_10,R0))' <-- This is another way of writing it with
% the lqr command

Ac = [ 0 0 0 0;
    L(:,1) (Aw - Bw*Kc - L*C)]

Bc1 = [1 0 0 0 0;
    zeros(3,2) L(:,2) zeros(3,2)]

Bc2 = [-1; -1; 0 ; 0]

Cc = [ 0 -Kc]

Dc1 = zeros(1,5)

Dc2 = 0

% To form a closed loop system, the controller is connected to the plant

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
     
disp('LQE; rho = 10^2 margins')
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

LQE2.y      = y;
LQE2.t      = t;
LQE2.az     = y(:,1); %  acceleration (fps2)
LQE2.q      = y(:,2); %  pitch rate (dps)
LQE2.dele   = y(:,3); %  dele (deg)
LQE2.dd     = y(:,4); %  dele_dot (dps)
LQE2.taur   = taur;
LQE2.taus   = taus;
LQE2.Lu     = Lu;
LQE2.rd     = rd;
LQE2.rdu_min   = rdu_min;
LQE2.sr     = sr;
LQE2.sru_min   = sru_min;
LQE2.S      = S;
LQE2.T      = T;
LQE2.lgcf   = wc;
LQE2.dele_Az  = dele_Az;
LQE2.dele_q   = dele_q;
LQE2.deledot_Az  = deledot_Az;
LQE2.deledot_q   = deledot_q;

step(sys_cl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Az comparison
figure('Name','Acceleration Time History')
plot(RSLQR.t,RSLQR.az,'LineWidth',3);
hold on
plot(LQE10.t,LQE10.az,'LineWidth',3);
plot(LQE4.t,LQE4.az,'LineWidth',3);
plot(LQE2.t,LQE2.az,'LineWidth',3);
xlabel('Time (sec)');
ylabel('Az (fps2)');
legend(['RSLQR'],...
    ['LQE; rho = 10^10'],...
    ['LQE; rho = 10^4'],...
    ['LQE; rho = 10^2'],...
    'Location','Best');
hold off


% Pitch Rate (q) comparison
figure('Name','Pitch Rate Time History')
plot(RSLQR.t,RSLQR.q,'LineWidth',3);
hold on
plot(LQE10.t,LQE10.q,'LineWidth',3);
plot(LQE4.t,LQE4.q,'LineWidth',3);
plot(LQE2.t,LQE2.q,'LineWidth',3);
xlabel('Time (sec)');
ylabel('Pitch Rate q (dps)');
legend(['RSLQR'],...
    ['LQE; rho = 10^10'],...
    ['LQE; rho = 10^4'],...
    ['LQE; rho = 10^2'],...
    'Location','Best');
hold off

% Dele comparison
figure('Name','Dele Time History')
plot(RSLQR.t,RSLQR.dele,'LineWidth',3);
hold on
plot(LQE10.t,LQE10.dele,'LineWidth',3);
plot(LQE4.t,LQE4.dele,'LineWidth',3);
plot(LQE2.t,LQE2.dele,'LineWidth',3);
xlabel('Time (sec)');
ylabel('Dele (deg)');
legend(['RSLQR'],...
    ['LQE; rho = 10^10'],...
    ['LQE; rho = 10^4'],...
    ['LQE; rho = 10^2'],...
    'Location','Best');
hold off

% Dele_dot comparison
figure('Name','Dele_dot Time History')
plot(RSLQR.t,RSLQR.dd,'LineWidth',3);
hold on
plot(LQE10.t,LQE10.dd,'LineWidth',3);
plot(LQE4.t,LQE4.dd,'LineWidth',3);
plot(LQE2.t,LQE2.dd,'LineWidth',3);
xlabel('Time (sec)');
ylabel('Deledot (dps)');
legend(['RSLQR'],...
    ['LQE; rho = 10^10'],...
    ['LQE; rho = 10^4'],...
    ['LQE; rho = 10^2'],...
    'Location','Best');
hold off

% Bode Mag Comparison
figure('Name','Bode Magnitude at Plant Input'),
semilogx(w,20*log10(abs(RSLQR.Lu)),...
    w,20*log10(abs(LQE10.Lu)),...
    w,20*log10(abs(LQE4.Lu)),...
    w,20*log10(abs(LQE2.Lu)),'Linewidth',3);
legend(['RSLQR'],...
    ['LQE; rho = 10^10'],...
    ['LQE; rho = 10^4'],...
    ['LQE; rho = 10^2'],...
    'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag')
title('Bode at Input')

% Plot Nyquist of Lu
% Need to define a unit circle for nyquist plot
n = 36; %// Define number of points on circle
theta = linspace(0, 2*pi, n);
x = cos(theta);
y = sin(theta);
figure('Name','Nyquist Plot at Plant Input'),
plot(x-1, y, '-.r'); % Unit circle
hold on
plot(real(squeeze(RSLQR.Lu)),imag(squeeze(RSLQR.Lu)),'Linewidth',2);
plot(real(squeeze(LQE10.Lu)),imag(squeeze(LQE10.Lu)),'Linewidth',2);
plot(real(squeeze(LQE4.Lu)),imag(squeeze(LQE4.Lu)),'Linewidth',2);
plot(real(squeeze(LQE2.Lu)),imag(squeeze(LQE2.Lu)),'Linewidth',2);
hold off
axis([-2 2 -2 2]);
legend('Unit Circle','RSLQR','LQE; rho = 10^10','LQE; rho = 10^4','LQE; rho = 10^2','Location','Best');
xlabel('Re(L)')
ylabel('Im(L)')
title('Nyquist Plot at Plant Input')


figure('Name','Return Difference at Plant Input'),
semilogx(w,20*log10(abs(RSLQR.rd)),...
    w,20*log10(abs(LQE10.rd)),...
    w,20*log10(abs(LQE4.rd)),...
    w,20*log10(abs(LQE2.rd)),'Linewidth',2)
legend(['RSLQR'],...
    ['LQE; rho = 10^10'],...
    ['LQE; rho = 10^4'],...
    ['LQE; rho = 10^2'],...
    'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Return Difference at Plant Input')


figure('Name','Stability Robustness at Plant Input'),
semilogx(w,20*log10(abs(RSLQR.sr)),...
    w,20*log10(abs(LQE10.sr)),...
    w,20*log10(abs(LQE4.sr)),...
    w,20*log10(abs(LQE2.sr)),'Linewidth',2)
legend(['RSLQR'],...
    ['LQE; rho = 10^10'],...
    ['LQE; rho = 10^4'],...
    ['LQE; rho = 10^2'],...
    'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Stability Robustness at Plant Input')

% Complementary Sensitivity plot
figure('Name','Comp Sensitivity T=Y/R')
semilogx(w,20*log10(abs(squeeze(RSLQR.T(1,1,:)))),'LineWidth',2);grid
hold on
semilogx(w,20*log10(abs(squeeze(LQE10.T(1,1,:)))),'LineWidth',2);grid

semilogx(w,20*log10(abs(squeeze(LQE4.T(1,1,:)))),'LineWidth',2);grid

semilogx(w,20*log10(abs(squeeze(LQE2.T(1,1,:)))),'LineWidth',2);grid

legend(['RSLQR'],...
    ['LQE; rho = 10^10'],...
    ['LQE; rho = 10^4'],...
    ['LQE; rho = 10^2'],...
    'Location','Best');
xlabel('Frequency (rps)')
ylabel('Co-Sens Mag (dB)')
title('Comp Sensitivity T=Y/R')
hold off

% Sensitivity plot
figure('Name',' Sensitivity S=E/R')
semilogx(w,20*log10(abs(squeeze(RSLQR.S(1,1,:)))),'LineWidth',2);grid
hold on
semilogx(w,20*log10(abs(squeeze(LQE10.S(1,1,:)))),'LineWidth',2);grid

semilogx(w,20*log10(abs(squeeze(LQE4.S(1,1,:)))),'LineWidth',2);grid

semilogx(w,20*log10(abs(squeeze(LQE2.S(1,1,:)))),'LineWidth',2);grid

legend(['RSLQR'],...
    ['LQE; rho = 10^10'],...
    ['LQE; rho = 10^4'],...
    ['LQE; rho = 10^2'],...
    'Location','Best');
xlabel('Frequency (rps)')
ylabel('Co-Sens Mag (dB)')
title('Sensitivity S=E/R')
hold off


% Plot Noise-To-Control
figure('Name','Noise-2-Control');
semilogx(w,RSLQR.dele_Az,'k',w,RSLQR.dele_q,'k--',...
         w,LQE10.dele_Az,'b',w,LQE10.dele_q,'b--',...
         w,LQE4.dele_Az,'g',w,LQE4.dele_q,'g--',...
         w,LQE2.dele_Az,'m',w,LQE2.dele_q,'m--',...
         'LineWidth',2);
legend('Az RSLQR','q RSLQR',...
       'Az LQE;rho=10^10','q LQE;rho=10^10',...
       'Az LQE;rho=10^4','q LQE;rho=10^4',...
       'Az LQE;rho=10^2','q LQE;rho=10^2',...
       'Location','Best');
title(['Noise to Control ']);
xlabel('Freq (rps)');ylabel('Mag (dB)');
grid;


% Plot Noise-To-Control Rate
figure('Name','Noise-2-Control Rate');
semilogx(w,RSLQR.deledot_Az,'k',w,RSLQR.deledot_q,'k--',...
         w,LQE10.deledot_Az,'b',w,LQE10.deledot_q,'b--',...
         w,LQE4.deledot_Az,'g',w,LQE4.deledot_q,'g--',...
         w,LQE2.deledot_Az,'m',w,LQE2.deledot_q,'m--',...
         'LineWidth',2);
legend('Az RSLQR','q RSLQR',...
       'Az LQE;rho=10^10','q LQE;rho=10^10',...
       'Az LQE;rho=10^4','q LQE;rho=10^4',...
       'Az LQE;rho=10^2','q LQE;rho=10^2',...
       'Location','Best');
title(['Noise to Control Rate']);
xlabel('Freq (rps)');ylabel('Mag (dB)');
grid;