%
% Adaptive OBLTR Design and Analysis
% K.A. Wise
% Created : 4/4/2012, Kevin A Wise
%
% Modified:
% 4/12/2017 K.A. Wise.   Adaptive
% LQR Q-Matrix line 113
% Observer vee and Qv Matrix line 216,221
% Adaptive Scaling, Deadzone line 446
% Uncertainty Modeling line 415
% y_cmd line 550
% System dynamics line 737
%%
clear all
close all
clc

format short e

disp('****************************** Program Start ****************');
plot_file_name = 'Hmk9_A_OBLTR_2017.ppt';
save_plots = 0; % Flag to bypass saving

rtd  = 180/pi;
d2r  = 1/rtd;
% g    = 32.174;
% dt = 0.002;
% t = 0:dt:2;
% w = logspace(-2,3,500);
% dd=0.:.001:2*pi;
% xx1=cos(dd)-1;yy1=sin(dd);

%*******************************************
% Define the Plant Model Ap and Bp
% Will design using second order model
% with AOA and q as states, dele the control
% Outputs are Az and q.
%*******************************************
% Public release airframe parameters

Za_V = -1.3046;
Ma   = 47.711; % Positve Ma = unstable aircraft
Zd_V = -.2142;
Md   = -104.83;
Mq   = 0.;
V    = 886.78; % (fps)
Za   = V*Za_V;
Zd   = V*Zd_V;
w_act = 2.*pi*50; % actuator natural frequency rps
z_act = 0.707;    % actuator damping
grav = 32.174; % fps2

disp('Actatuor wn Hz')
disp(w_act/2/pi)

% Define the nominal plant matrices
% **************************************************************************
% Plant model for analysis with actuator 
% **************************************************************************
% States are AOA (rad),  pitch rate (rps), dele (rad), deledot (rps)
Ap = [Za_V  1.      Zd_V          0.; 
      Ma    Mq       Md           0.;
      0.    0.       0.           1.;
      0.    0. -w_act*w_act -2*z_act*w_act];
% Input is dele (rad)
Bp = [0.; 0.; 0.; w_act*w_act ];

% Outputs are Az (fps2), AOA (rad), pitch rate (rps), dele (rad), deledot (rps)
Cp = [  Za   0.  Zd  0.; 
        eye(4)];
Dp = 0.*Cp*Bp; 

% State space model for the actuator
A_act = [  0.           1.;
       -w_act*w_act -2*z_act*w_act];
B_act = [ 0.; w_act*w_act];
C_act = [ 1 0];
D_act = 0.;
sys_act = ss(A_act,B_act,C_act,D_act);
  

% plant size
[nCp,nAp] = size(Cp);
[~,nBp]   = size(Bp);
mp = 1;
mBp = 1;
nCp_reg = 1;
sel_state_alpha_sp = 1;
sel_state_q_sp     = 2;


%%
%*******************************************
% Design the RSLQR using 3rd order Aw and Bw
% No actuator in the design model
%*******************************************
% State Feedback Design
% states are [ int-err AOA q ]
% zdot = Aw*z + Bw*mu + Bcmd*r

Aw = [ 0.   Za  0. ;
       0.   Za_V  1. ;
       0.   Ma    Mq ];
Bw = [ Zd;
       Zd_V;
       Md];
[nw, mw] = size(Bw);   
Bcmd = [ -1.;
          0.;
          0.];  
    
Cw = [ 0.  Za  0.];  % Form Az from the wiggle system
Dw = Zd;             % Form Az from the wiggle system

%--------------------------------------
% Insert your LQR Q matrix here
%--------------------------------------
%LQR penalty matrix Q
Q_LQR = 0.*eye(3);
Q_LQR(1,1) = 0.9374e-04;
R_LQR = 1.;

% Solve for RSLQR state feedback gains
[Kx_lqr,Pw]=lqr(Aw,Bw,Q_LQR,R_LQR);

Kw = Kx_lqr;

Aw_cl = Aw - Bw*Kw;
Cw_cl = Cw - Dw*Kw;
Dw_cl = Dw*0;


% LQR solution check
lqr_eq = Pw*Aw + Aw'*Pw + Q_LQR - Pw*Bw*inv(R_LQR)*Bw'*Pw;
norm_lqr_eq = norm(lqr_eq);
disp(['LQR ARE Solution Tolerance is: ', num2str(norm_lqr_eq)]);

%*************************************************************************
% Analyze LQR Design
%*************************************************************************

% populate the controller matrices  
% y = Az (fps2), AOA (rad), pitch rate (rps), dele (rad), deledot (rps)
% u = - Kx_lqr*[ int-err aoa, q ];
    Ac =  0.;
    Bc1 = [1. 0. 0. 0. 0.];
    Bc2 =  -1;
    Cc = -Kx_lqr(1);
    Dc1 = [0. -Kx_lqr(2:3) 0. 0.];
    Dc2 = 0.;

% Form the closed loop system
     Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
     Acl = [     (Ap+Bp*Z*Dc1*Cp)              (Bp*Z*Cc);
         (Bc1*(Cp+Dp*Z*Dc1*Cp))  (Ac+Bc1*Dp*Z*Cc)];
     Bcl = [       Bp*Z*Dc2;
         (Bc2+Bc1*Dp*Z*Dc2)];
     Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
     Dcl =(Dp*Z*Dc2);
     sys_cl = ss(Acl,Bcl,Ccl,Dcl);

 % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 % Closed-loop stability check
 % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     if max(real(eig(Acl))) > 0,
         disp('Closed-Loop System is Unstable');
         disp(['  Most unstable eig = ', num2str(max(real(eig(Acl))))]);
         return
     end;          
      
 disp('Homework 3 LQR Design')
 disp('[Aw Bw]')
 disp([Aw Bw])
 disp('Q_LQR')
 disp(Q_LQR)
 disp('R_LQR')
 disp(R_LQR)
 disp('Kx_lqr')
 disp(Kx_lqr)
 disp('  ')
   
     
     
     
%%
%*************************************************************************
%   OBLTR
%*************************************************************************
%  observer states = {'int-err' 'aoa' 'q'};
%  control = {'dele' };
%  measurements = {'int-err' 'q' };

Aw_lqg = Aw;
Bw_lqg = Bw;
Cw_lqg = [ 1 0 0;  % int-error
           0 0 1]; % pitch rate
Dw_lqg = [ 0; 0];

Cw_meas = Cw_lqg;
Dw_meas = Dw_lqg;

Cw_all_cl = Cw_meas - Dw_meas*Kw;
Dw_all_cl = Dw_meas*0;

Cw_meas_cl = Cw_meas - Dw_meas*Kw;
Dw_meas_cl = Dw_meas*0;

[nCw_meas, mCw_meas] = size(Cw_meas);


% Plant disturbance and measurement noise covariance matrices
Q0 = diag([ 0.001 0.0014 0.005 ]);
R0 = 1000*diag([ 0.025^2 0.001^2]);

zz =  -10 ;% Place this zero when squaring up the system

q_scale = 0.5;

%  Square Up The Dynamics
sys_1 = ss(Aw_lqg, Bw_lqg, Cw_lqg, Dw_lqg);
sys_2 = square_system_inputs_pp(sys_1,zz);
disp(['Bbar Matrix 1 z = ' num2str(zz) ])
sys_2.b    
tzero(sys_2)

nQ = norm(Q0);
nBbar = norm(sys_2.b*sys_2.b');
[~,nAw] = size(Aw);


%--------------------------------------
% Insert your LTR vee  here
%--------------------------------------

vee =  .0004;

%--------------------------------------
% Insert your Qv,Rv for observer gain design here
%--------------------------------------
Qv = Q0 + q_scale*(nQ/nBbar)*((vee+1)/vee)*(sys_2.b*sys_2.b');
Rv = (vee/(vee+1))*R0;

[Kf_lqg,Pv] = lqe(Aw_lqg,eye(3),Cw_lqg,Qv,Rv);


[U, Lambda, V] = svd(sys_2.b'*sys_2.c'*inv(sqrt(R0)));
W = (U*V')';
% PB = C testing
PwBw = inv(Pv)*sys_1.b;
Cw_computed = sys_2.c'*inv(sqrt(R0))*W;
    Cw_computed = Cw_computed(:,1)';
    Check_PBC = [PwBw, Cw_computed'];
    disp('PwBwCw Data');
    disp(Check_PBC);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Directions of adaptations, (need to be small)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sqrt_diag_Cw_computed_squared = sqrt(diag(Cw_computed*Cw_computed'));
sqrt_diag_PwBw_squared        = sqrt(diag(PwBw'*PwBw));
diag_Cw_computedPwBw          = diag((Cw_computed*PwBw));
angle_PBC_deg = acosd(diag_Cw_computedPwBw./sqrt_diag_Cw_computed_squared./sqrt_diag_PwBw_squared);
% angle_PBC_deg = acosd(diag(Cw_computed*PwBw/norm(Cw_computed)/norm(PwBw)));
disp(  'Adaptation Directions Angle (deg)');
disp(num2str(angle_PBC_deg));

Kwf = Kf_lqg';

% KF solution check - important to check this when using LTR
kf_eq = Aw_lqg*Pv + Pv*Aw_lqg' + Qv - Pv*Cw_lqg'*inv(Rv)*Cw_lqg*Pv;
norm_kf_eq = norm(kf_eq);
disp(['KF FARE Solution Tolerance is: ', num2str(norm_kf_eq)]);

% LQG Controller Model
Ac = [       0.*ones(1,4)      ;
    Kf_lqg(:,1)  (Aw - Bw*Kx_lqr - Kf_lqg*Cw_lqg)];
Bc1 = [   1.         0.  0.         0.  0.  ;
         0.*ones(nAw,2) Kf_lqg(:,2) 0.*ones(nAw,2) ];
Bc2 = [ -1;
        -1;
    0.*ones(2,1)];
Cc = [0. -Kx_lqr];
Dc1 = [ 0. 0. 0. 0. 0.];
Dc2 = 0.;
    
% Form the closed loop system
     Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
     Acl = [     (Ap+Bp*Z*Dc1*Cp)              (Bp*Z*Cc);
         (Bc1*(Cp+Dp*Z*Dc1*Cp))  (Ac+Bc1*Dp*Z*Cc)];
     Bcl = [       Bp*Z*Dc2;
         (Bc2+Bc1*Dp*Z*Dc2)];
     Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
     Dcl =(Dp*Z*Dc2);
     sys_cl = ss(Acl,Bcl,Ccl,Dcl);

 % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 % Closed-loop stability check
 % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     if max(real(eig(Acl))) > 0,
         disp('Closed-Loop System is Unstable');
         disp(['  Most unstable eig = ', num2str(max(real(eig(Acl))))]);
         return
     end;     

disp('Observer Design')
 disp('[Aw Cw^T]')
 disp([Aw_lqg Cw_lqg'])
 disp('vee')
 disp(vee)
 disp('Qv')
 disp(Qv)
 disp('Rv')
 disp(Rv)
 disp('Lv')
 disp(Kf_lqg)


     
%%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Adaptive Control Simulation 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Redfine Plant to just have short period Dynamics
Ap = Ap( 1:2, 1:2);
Bp = Bp( 1:2);
Cp = eye(2);
Dp = 0.*Cp*Bp;


uncertainty_choice = menu('Plant Uncertainty:','Off','Subtract u_bl', 'Subtract u_bl + RBF');
control_choice     = menu('Control Choice:','LQR','OBLTR','OBLTR + Adaptive');
actuator_choice    = menu('Control Actuators:','Off','On') - 1;

disp(' ');
if uncertainty_choice == 1,
  disp('Plant Uncertainty Off');
  U_title = 'Unc = OFF, ';
elseif uncertainty_choice == 2,
  disp('Linear Plant Uncertainty selected');
  U_title = 'Unc = ubl, ';
  elseif uncertainty_choice == 3,
  disp('(Linear + RBF) Plant Uncertainty selected');
  U_title = 'Unc = ubl+RBF, ';
end;

if control_choice == 1,
  disp('LQR State Feedback selected');
  
  lqr_flag            = 1;  % LQR is On
  obltr_flag          = 0;  % OBLTR is Off
  obltr_adaptive_flag = 0;  % OBLTR + Adaptive is Off
  C_title = 'K = LQR, ';

elseif control_choice == 2,
  disp('OBLTR Output Feedback selected');
  
  lqr_flag            = 0;  % LQR is Off
  obltr_flag          = 1;  % OBLTR is On
  obltr_adaptive_flag = 0;  % OBLTR + Adaptive is Off
  C_title = 'K = OBLTR, ';

elseif control_choice == 3,
  disp('(OBLTR + Adaptive) Output Feedback selected');
  
  lqr_flag            = 0;  % LQR is Off
  obltr_flag          = 0;  % OBLTR is Off
  obltr_adaptive_flag = 1;  % OBLTR + Adaptive is On
  C_title = 'K = OBLTR+Adaptive, ';

end;

if actuator_choice == 0,
  disp('Actuators Off');
  A_title = 'Act = OFF ';

else
  disp('Actuators On');
  A_title = 'Act = ON ';
end;

plot_title = [ U_title C_title A_title];

drawnow

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Simulation time data
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dt     = 0.001;  % integration time step, sec
Tfinal = 50;    % simulation stop time , sec

time   = 0:dt:Tfinal;
npnts  = length(time);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Discretization data, (for simulation only)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sysw_ol  = ss(Aw, [Bw, Bcmd, Kwf'], Cw, [Dw, zeros(mp,mp+nCw_meas)]);
[syswd_ol, T_c2d_ol] = c2d(sysw_ol,dt);

Aw_d   = syswd_ol.a;
Bw_d   = syswd_ol.b(:,1:mp);
Bcmd_d = syswd_ol.b(:,mp+1:mp+mp);
Kwf_d  = syswd_ol.b(:,mp+mp+1:end)';

Cw_d   = syswd_ol.c;
Dw_d   = syswd_ol.d;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Reference model data for adaptive control
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Aref = Aw_cl;
Bref = Bcmd;
Cref = Cw_cl;
Dref = Dw_cl;

sysw_ref  = ss(Aref, Bref, Cref, Dref);
[syswd_ref, T_c2d_ref] = c2d(sysw_ref,dt);

Aref_d   = syswd_ref.a;
Bref_d   = syswd_ref.b;
Cref_d   = syswd_ref.c;
Dref_d   = syswd_ref.d;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Actuator
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[sys_act_d, T_c2d_dele] = c2d(sys_act, dt);

A_act_d = sys_act_d.a;
B_act_d = sys_act_d.b;
C_act_d = sys_act_d.c;
D_act_d = sys_act_d.d;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Discrete closed-loop system stability check
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Aw_d_cl = Aw_d - Bw_d*Kw;  % Discrete closed-loop matrix

if max(abs(eig(Aw_d_cl))) > 1,
  disp('Discrete Closed-Loop System is Unstable');
  disp(['  Most unstable eig = ', num2str(max(real(eig(Aw_d_cl))))]);
  disp('Simulation Stopped');
  return
else
  disp('Discrete Closed-Loop System is Stable');
  disp(['  System Dominant Eig = ', num2str(max(real(eig(Aw_d_cl))))]);
  disp(['  System Fast     Eig = ', num2str(min(real(eig(Aw_d_cl))))]);
end;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% System uncertainty
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Lambda = eye(mp);
dKw = Kw*0;

if uncertainty_choice > 1,
  Lambda      = 1.0;                     % loss in control effectiveness
  %dKw(:,mp+1:end) = Kw(:,mp+1:end)*1.0;  % incremental loss in proportional feedback gains
  dKw = Kw;  % Subtracts u_baseline
end;

% Nonlinear-in-alpha uncertainty
alpha_data      = [-4:.1:4]*d2r;
alpha_data      = [-6:.1:6]*d2r;
alpha_data_min  = min(alpha_data);
alpha_data_max  = max(alpha_data);

RBF_sigma_alpha = 2/3*1.0*d2r;          % RBF sigma
RBF_W_alpha     = 1/RBF_sigma_alpha^2;  % RBF weight
RBF_center      = 2.0*pi/180       ;    % RBF center

d_Fxw_alpha_data = rbf2mat_v2(alpha_data,RBF_center,RBF_W_alpha*2)*0.3;
d_Fxw_alpha_data = d_Fxw_alpha_data - d_Fxw_alpha_data(1);

if uncertainty_choice ~= 3,
  d_Fxw_alpha_data = d_Fxw_alpha_data*0;
end;

Fx = 0;  % modeling uncertainty IC

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Adaptation rate scale  = O(norm(Kwf))
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%AdaptiveRate_scale = 10.0*norm(Kwf);
AdaptiveRate_scale = 0.008*norm(Kwf);  % Scales the learning rates
dead_zone          = 1e-06;
%dead_zone = 0.04;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% KX adaptation rates for linear regressor [u_bl; xw_hat]. Set to norm(Kwf)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Gamma_u_bl = AdaptiveRate_scale*eye(mBp);  % Adaptation rates for baseline control

% Adaptation rates  for observer states
Gamma_eyI       = zeros(mBp);
Gamma_alpha_hat = 1;
Gamma_q_hat     = 1;

% Total adaptation rate matrix for linear regressor
Gamma_X = zeros(mBp+nw);

Gamma_X(1:mBp,1:mBp)             = Gamma_u_bl;
Gamma_X(mBp+1:2*mBp,mBp+1:2*mBp) = Gamma_eyI;
Gamma_X(2*mBp+1:end,2*mBp+1:end) = diag([Gamma_alpha_hat, Gamma_q_hat])*AdaptiveRate_scale;

gamma_X = 0;   % e-mod gain for Kxw_hat
Kxw_max = 10*ones(1,mBp);  % projection bounds for Kxw_hat, columns-wise

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% RBF break points
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha_bp_min = -2*d2r;    % minimum break point,   rad
alpha_bp_max =  2*d2r;    % maximum break point,   rad
d_alpha_bp   =  1/2*d2r;  % break point increment, rad

alpha_bp   = [alpha_bp_min:d_alpha_bp:alpha_bp_max]; % array of break points, rad
nalpha_bp  = length(alpha_bp);  % number of break points

N_rbf = nalpha_bp;   % total number of RBFs

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% RBF centers: RBF = function(alpha)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RBF_centers = zeros(nw,N_rbf);
k = 0;
for iv = 1:nalpha_bp,
  k = k + 1;
  RBF_centers(mBp+1,k) = alpha_bp(iv);
end;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% RBF sigma-weights
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RBF_W = zeros(nw,nw);  % IC matrix of RBF width

RBF_sigma_alpha = 4/3*d_alpha_bp;  % RBF sigma for alpha_bp
RBF_W(mBp + sel_state_alpha_sp,mBp + sel_state_alpha_sp) = 1/RBF_sigma_alpha^2;  % RBF weight

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot RBF = RBF(alpha) data
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha_range_data = [-4:0.1:4]*d2r;
X_range = zeros(nw, length(alpha_range_data));
X_range(mBp + sel_state_alpha_sp,:) = alpha_range_data;

RBF_range_data = rbf2mat_v2(X_range, RBF_centers, RBF_W);

norm_RBF_range_data = zeros(length(alpha_range_data),1);
for kk = 1:length(alpha_range_data);
  norm_RBF_range_data(kk) = norm(RBF_range_data(kk,:));
end;

figure('Name','RBF Model'),
subplot(211);
plot(alpha_range_data*rtd, RBF_range_data); grid;
title('RBF Data');
ylabel('RBF(alpha)');

subplot(212);
plot(alpha_range_data*rtd, norm_RBF_range_data); grid;
ylabel('RBF(alpha)');
xlabel('alpha, deg');
if(save_plots == 1) saveppt2(plot_file_name); end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% RBF adaptation rates = norm(Kwf)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Gamma_Theta  = AdaptiveRate_scale*eye(N_rbf+1);  % Theta_hat learning rates
Gamma_Theta(end,end) = 0;        % Bias adaptation rate for Theta_hat is Off

gamma_Theta = 0;   % Theta_hat e-mod gain

Theta_max  = 10.0*ones(1,mBp); % projection bounds for Theta_hat, columns-wise

N = N_rbf + 1;     % regressor dimension, includes bias

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Other data
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MixY = eye(mp,mp);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% External Commands
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
y_cmd_max = -100;   % Az cmd , fps2

% Series of doublets
t_cmd_table = [0.0, 1.0, 1.1, 10.0, 10.1, 15.0, 15.1,  25.0, 25.1, 30.0, 30.1, 40.0, 40.1, Tfinal]; 
y_cmd_table = [0.0, 0.0  1.0,  1.0,  0.0,  0.0,  1.0,  1.0,  0.0,   0.0,  1.0,  1.0,  0.0,  0.0  ];

% Command time constant
s = tf('s');
tau_cmd = 0.2;  % 200 Msec
sys_cmd = ss(1/(tau_cmd*s+1));

y_cmd_interp = interp1(t_cmd_table,y_cmd_table,time);

% y_cmd_sim = lsim(sys_cmd,y_cmd_interp,time)*y_cmd_max;
% y_cmd = y_cmd_sim';
y_cmd = y_cmd_max*y_cmd_interp;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Array data pre-allocation
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
xw = zeros(nw,npnts);
% xw(mp+1:nw,1) = 0.02*(-1 + 2*rand(nw-mp,1));

xw_hat = zeros(nw,npnts);
% xw_hat(mp+1:nw,1) = 0.02*(-1 + 2*rand(nw-mBp,1));

xw_ref = zeros(nw,npnts);
% xw_hat(:,1) = -1 + 2*rand(nw,1);  % testing observer performance

xw_dot     = zeros(nw,npnts);
xw_hat_dot = zeros(nw,npnts);
xw_ref_dot = zeros(nw,npnts);

norm_e_data   = zeros(npnts,1);
norm_exw_data = zeros(npnts,1);
norm_ey_data  = zeros(npnts,1);

u_cmd = zeros(mBp,npnts);
u     = zeros(mBp,npnts);
u_bl  = zeros(mBp,npnts);
u_ad  = zeros(mBp,npnts);

% RBF parameters
Theta_hat           = zeros(N,mBp); 
Theta_hat_new       = zeros(N,mBp); 
Theta_hat_dot       = zeros(N,mBp); 
Theta_hat_dot_p     = zeros(N,mBp);
norm_Theta_hat_data = zeros(npnts,mBp);
Theta_hat_data      = zeros(npnts,N,mBp);

% Linear regressor [u_bl; xw_hat] parameters
Kxw_hat            = zeros(nw+mBp,mBp); 
Kxw_hat_new        = zeros(nw+mBp,mBp); 
Kxw_hat_dot        = zeros(nw+mBp,mBp); 
Kxw_hat_dot_p      = zeros(nw+mBp,mBp);

Kxw_hat_data       = zeros(npnts,nw+mBp,mBp);
norm_Kxw_hat_data  = zeros(npnts,1);

Kxw_hat_u_bl       = zeros(npnts,mBp,mBp);
norm_Kxw_hat_u_bl  = zeros(npnts,1);

Kxw_hat_xw_hat      = zeros(npnts,nw,mBp);
norm_Kxw_hat_xw_hat = zeros(npnts,1);

Fxw_data = zeros(mBp,npnts);

x_act = zeros(2*mBp,npnts);
y_act = zeros(mBp,npnts);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Simulation loop starts here
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp(' Simulation in progress ...');
drawnow

for k = 1:npnts-1,
  xw_k     = xw(:,k);      % system state
  xw_hat_k = xw_hat(:,k);  % observer state
  xw_ref_k = xw_ref(:,k);  % ref model state
  
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % Matched uncertainty
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  alpha_k = xw_k(mp+sel_state_alpha_sp);  % alpha, rad
  alpha_k = max( alpha_data_min, min(alpha_k, alpha_data_max) );
  dFxw_k = interp1(alpha_data,d_Fxw_alpha_data,alpha_k);
  
  % Total uncertainty
  Fxw_k         = dKw*xw_k + dFxw_k;
  Fxw_data(:,k) = Fxw_k;

  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % Baseline control
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if lqr_flag == 1,        
    u_bl_k = -Kw*xw_k;      % LQR state feedback
  else                     
    u_bl_k = -Kw*xw_hat_k;  % OBLTR feedback
  end;
  
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % Adaptive augmentation
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if obltr_adaptive_flag == 1,
    Y_k     = Cw_computed*xw_k ;    % system   output for adaptive laws
    Y_hat_k = Cw_computed*xw_hat_k; % observed output for adaptive laws

    [Kxw_hat_dot, Theta_hat_dot, ...
     Kxw_hat_new, Theta_hat_new, u_ad_k] = ...
        Proj_OFMRAC(Y_k, Y_hat_k, xw_hat_k, u_bl_k, ...
                    Kxw_hat, Theta_hat, ...
                    Kxw_hat_dot_p, Theta_hat_dot_p,  ...
                    RBF_centers, RBF_W, MixY, ...
                    Gamma_X, Gamma_Theta, ...
                    gamma_X, gamma_Theta, dead_zone, ...
                    Kxw_max, Theta_max, dt);
     
    % current Kxw_hat and Theta_hat values
    Kxw_hat_data(k,:,:)    = Kxw_hat;
    Theta_hat_data(k,:,:)  = Theta_hat;
    norm_Kxw_hat_data(k)   = norm(Kxw_hat);
    norm_Theta_hat_data(k) = norm(Theta_hat);
    
    Kxw_hat_u_bl(k,:,:)    = Kxw_hat(1:mBp,:);
    norm_Kxw_hat_u_bl(k)   = norm(Kxw_hat(1:mBp,:));
    
    Kxw_hat_xw_hat(k,:,:)  = Kxw_hat(mBp+1:end,:);
    norm_Kxw_hat_xw_hat(k) = norm(Kxw_hat(mBp+1:end,:));
    
    % Kxw_hat and Theta_hat values for next iteration
    Kxw_hat       = Kxw_hat_new;
    Kxw_hat_dot_p = Kxw_hat_dot;
    
    Theta_hat       = Theta_hat_new;
    Theta_hat_dot_p = Theta_hat_dot;
  else
    u_ad_k = 0;
  end;
  
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % Control commands
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  u_cmd_k    = u_bl_k + u_ad_k;  % total control
  u_cmd(:,k) = u_cmd_k;
  
  % Total control input
  u(:,k) = u_cmd_k;
  u_k    = u(:,k);
  
  u_bl(:,k) = u_bl_k;   % baseline control
  u_ad(:,k) = u_ad_k;   % adaptive increment
  
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % Actutator dynamics propagation
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if actuator_choice == 1,
    x_act(:,k+1) = A_act_d*x_act(:,k) + B_act_d*u_cmd_k;
    y_act(:,k)   = C_act_d*x_act(:,k) + D_act_d*u_cmd_k;
    
    u(:,k) = y_act(:,k);
    u_k    = u(:,k);
  end;
  
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % Outputs
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  y_k = Cw_meas*xw_k + Dw_meas*Lambda*(u(:,k) + Fxw_k);  % system
  
  y_hat_k = Cw_meas_cl*xw_hat_k;  % observer
  y_ref_k = Cw_meas_cl*xw_ref_k;  % ref model

  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % Errors
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  e_k   = xw_hat_k - xw_ref_k;  % observer state  tracking   error
  exw_k = xw_hat_k - xw_k;      % observer state  estimation error
  ey_k  = y_hat_k - y_k;        % observer output estimation error
  
  norm_e_data(k)   = norm(e_k);
  norm_exw_data(k) = norm(exw_k);
  norm_ey_data(k)  = norm(ey_k);
  
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % Dynamics propagation
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % System
  xw(:,k+1) = Aw_d*xw_k + Bw_d*Lambda*(u(:,k) + Fxw_k) + Bcmd_d*y_cmd(:,k);
  
  % Observer
  xw_hat(:,k+1) = Aw_d_cl*xw_hat(:,k) + Bcmd_d*y_cmd(:,k) + Kwf_d'*(y_k - y_hat_k);

  % Ideal ref model
  xw_ref(:,k+1) = Aref_d*xw_ref(:,k) + Bref_d*y_cmd(:,k);
end;
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Main simulation loop ends here
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

u_cmd(:,end) = u_cmd(:,end-1);
u(:,end)     = u(:,end-1);
u_bl(:,end)  = u_bl(:,end-1);
u_ad(:,end)  = u_ad(:,end-1);

norm_Kxw_hat_data(end)   = norm_Kxw_hat_data(end-1);
norm_Theta_hat_data(end) = norm_Theta_hat_data(end-1);

norm_e_data(end)   = norm_e_data(end-1);
norm_exw_data(end) = norm_exw_data(end-1);
norm_ey_data(end)  = norm_ey_data(end-1);

Kxw_hat_data(end,:)   = Kxw_hat_data(end-1,:);
Kxw_hat_u_bl(end,:)   = Kxw_hat_u_bl(end-1,:);
Kxw_hat_xw_hat(end,:) = Kxw_hat_xw_hat(end-1,:);

norm_Kxw_hat_u_bl(end)  = norm_Kxw_hat_u_bl(end-1);
norm_Kxw_hat_xw_hat(end) = norm_Kxw_hat_xw_hat(end-1);

% Lambda hat calculation
Lambda_hat_dele = 1./(1 - Kxw_hat_u_bl);

% Control rates
[Y_act, T_act, X_act] = lsim(sys_act, u_cmd, time);
u_dot = X_act(:,2)';

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% States
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
xp     = xw(mw+1:end,:);      % system 
xp_hat = xw_hat(mw+1:end,:);  % estimated 
xref   = xw_ref(mw+1:end,:);  % ideal ref model 

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Outputs
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
y = Cw*xw + Dw*Lambda*(u + Fxw_data);  % system regulated output

y_hat = Cw_cl*xw_hat;  % observer  regulated output
y_ref = Cref*xw_ref;   % ref model regulated output

% System, observer, and ref model output
y_all     = Cp*xw(nCp_reg+1:end,:)     + Dp*Lambda*(u + Fxw_data);
y_hat_all = Cp*xw_hat(nCp_reg+1:end,:) + Dp*u_bl;

y_ref_all = [xw_ref(mp+sel_state_alpha_sp,:);
             xw_ref(mp+sel_state_q_sp,:);
             y_ref];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Controls
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dele_cmd_deg  = u_cmd(1,:)*rtd;  % dele cmd 
dele_deg      = u(1,:)*rtd;      % dele 
dele_dot_dps  = u_dot(1,:)*rtd;  % dele rate

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Output signals
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha_deg     = xw(nCp_reg+sel_state_alpha_sp,:)*rtd;
alpha_hat_deg = xw_hat(nCp_reg+sel_state_alpha_sp,:)*rtd;
alpha_ref_deg = xw_ref(nCp_reg+sel_state_alpha_sp,:)*rtd;

q_dps     = xw(nCp_reg+sel_state_q_sp,:)*rtd;
q_hat_dps = xw_hat(nCp_reg+sel_state_q_sp,:)*rtd;
q_ref_dps = xw_ref(nCp_reg+sel_state_q_sp,:)*rtd;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Commands
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%alpha_cmd_deg = y_cmd*rtd;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Command tracking
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure('Name','Command Following'),
plot(time, y_cmd, time, y_ref, time, y, time, y_hat,'-.', 'Linewidth',2), grid;
ylabel('Az, fps2','fontsize',10);
title(['Tracking Performance ' plot_title],'fontsize',10);
legend('Cmd','Reference','Actual','Estimated','fontsize',10,'Location','Best');
xlabel('Time, sec','fontsize',10);
if(save_plots == 1) saveppt2(plot_file_name); end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Control inputs and rates
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure('Name','Control Usage'),
subplot(211);
if actuator_choice == 0,
  plot(time, dele_deg, 'Linewidth',2),grid,
else
  plot(time, dele_cmd_deg, time, dele_deg, 'Linewidth',2),grid,
  legend('Cmd', 'Actual', 'fontsize',10);
end;
title(['Elevon Pos and Rate ' plot_title],'fontsize',10);
ylabel('Dele','fontsize',10);
subplot(212);
plot(time, dele_dot_dps,'Linewidth',2),grid,
ylabel('Dele dot','fontsize',10);
xlabel('Time, sec','fontsize',10);
if(save_plots == 1) saveppt2(plot_file_name); end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% States
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure('Name','States AOA Q'),
subplot(211);
plot(time, alpha_ref_deg, time, alpha_deg, time, alpha_hat_deg, 'Linewidth',2), grid;
title(['States ' plot_title],'fontsize',10);
ylabel('alpha, deg','fontsize',10);
legend('Reference', 'Actual', 'Estimated', 'fontsize',10);
subplot(212);
plot(time, q_ref_dps, time, q_dps, time, q_hat_dps, 'Linewidth',2), grid;
ylabel('q, dps','fontsize',10);
if(save_plots == 1) saveppt2(plot_file_name); end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Adaptive parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure('Name','Adaptive Parameters'),
subplot(311);
plot(time, Lambda_hat_dele,'Linewidth',2), grid;
title(['Adaptive Parameter Norms ' plot_title],'fontsize',10);
ylabel('Lambda dela','fontsize',10);

subplot(312);
plot(time, norm_Kxw_hat_data,'Linewidth',2), grid;
ylabel('K x hat','fontsize',10);

subplot(313);
plot(time, norm_Theta_hat_data,'Linewidth',2), grid;
ylabel('Theta hat','fontsize',10);
xlabel('Time, sec','fontsize',10);
if(save_plots == 1) saveppt2(plot_file_name); end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Estimation and Tracking errors
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure('Name','Tracking errors'),
subplot(211);
plot(time, norm_exw_data,'Linewidth',2), grid;
title(['Estimation and Tracking Errors ' plot_title],'fontsize',10);
ylabel('||xhat - x||','fontsize',10);
subplot(212);
plot(time, norm_e_data,'Linewidth',2), grid;
ylabel('||xhat - xref||','fontsize',10);
xlabel('Time, sec','fontsize',10);
if(save_plots == 1) saveppt2(plot_file_name); end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Matched uncertainty approximation
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Fxw_hat_data = -(eye(mBp) - inv(Lambda))*u_bl - u_ad;

figure('Name','Uncertainty Approx'),
subplot(211);
plot(time, Fxw_data, time, Fxw_hat_data); grid;
title(['Dele: Matched Uncertainty ' plot_title],'fontsize',10);
legend('Actual','Estimated','fontsize',10);
xlabel('Time, sec','fontsize',10);
subplot(212);
plot(alpha_deg, Fxw_data, alpha_deg, Fxw_hat_data); grid;
xlabel('alpha, deg','fontsize',10);
if(save_plots == 1) saveppt2(plot_file_name); end

disp('Done');

return
