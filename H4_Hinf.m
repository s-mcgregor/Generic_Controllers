clear all
close all

plot_file_name = 'x45a_hinf_design.ppt';
save_plots = 0; % Flag to bypass saving plots

rtd  = 180/pi;
g    = 32.174;

m2ft = 3.2808;    % meters to feet conversion
ft2m = 1/3.2808;  % feet to meters conversion

Za_V = -1.3046;
Ma   = 47.711;
Zd_V = -.2142;
Md   = -104.83;
V    = 886.78;
Za   = V*Za_V;
Zd   = V*Zd_V;
Mq = -1.03341;
wa = 2.*pi*11; 
za = 0.707; 

%% Plant Model
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plant model with Az state
%   states   = [Az, q, deltae, deltae_dot]
%   controls = [deltaec]
%   outputs  = [Az, q, deltae, deltae_dot]
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Plant.Ap = [Za_V  1.      Zd_V          0.;

            Ma    0.       Md           0.;

            0.    0.       0.           1.;

                 0.    0. -wa*wa -2*za*wa];

Plant.Bp = [0.; 0.; 0.; wa*wa ]; 

Plant.Cp = [1,0,0,0];

Plant.Dp = 0.*Plant.Cp*Plant.Bp;
disp(['Open-Loop System Eigenvalues'])
damp(eig(Plant.Ap))

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~  Hinf Optimal Control  ~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Loop Gain Cross-Over Frequency

wlgcf = 22.3; %Hz, selected from HW3

%wlgcf = 1.4; %Hz

%% Sensitivity Ws

Ws.Tau = 1/(2*pi*wlgcf);
Ws.K = 0.5/Ws.Tau;
%Ws.K = 0.8/Ws.Tau;
Ws.num = Ws.K * [Ws.Tau, 1];
Ws.den = [1, 0];
[Ws.A,Ws.B,Ws.C,Ws.D] = tf2ss(Ws.num,Ws.den);

Ws.sys = ss(Ws.A,Ws.B,Ws.C,Ws.D);
figure(1)
bodemag(Ws.sys)
grid on
hold on

%% Complementary Sensitivity Wt

Wt.TauN = 1/(2*pi*wlgcf);
Wt.TauD = 0.005;
Wt.K = 0.707;
Wt.num = Wt.K * [Wt.TauN, 1];
Wt.den = [Wt.TauD, 1];
[Wt.A,Wt.B,Wt.C,Wt.D] = tf2ss(Wt.num,Wt.den);

Wt.sys = ss(Wt.A,Wt.B,Wt.C,Wt.D);
bodemag(Wt.sys)
grid on
hold on

%% Control Activity

Wc = 0.1;
Cc = [0, 0, 0, 1];
Wc_sys = ss(0,0,0,Wc);

bodemag(Wc_sys)
xlim([0.1,100])
if(save_plots == 1) saveppt2(plot_file_name); end   
hold off

%% Hinf State-Feedback Matrix Setup

HinfSF.A = [     Plant.Ap, zeros(4,1), zeros(4,1);
            -1*Ws.B*Plant.Cp,    Ws.A,          0;
            Wt.B*Plant.Cp,          0,       Wt.A];

HinfSF.B = [        Plant.Bp;
            -1*Ws.B*Plant.Dp;
               Wt.B*Plant.Dp];

HinfSF.E = [zeros(4,1);
                  Ws.B;
                     0];

HinfSF.C = [ -1*Ws.D*Plant.Cp, Ws.C,    0;
             Wt.D*Plant.Cp,    0, Wt.C;
            Wc*Cc*Plant.Ap,    0,    0];
            
HinfSF.D1 = [ -1*Ws.D*Plant.Dp;
              Wt.D*Plant.Dp;
             Wc*Cc*Plant.Bp];

HinfSF.D2 = [Ws.D;
                0;
                0];

HinfSF.B_hat = [HinfSF.B, HinfSF.E];

HinfSF.D_hat = [HinfSF.D1, HinfSF.D2];

HinfSF.Q = HinfSF.C'*HinfSF.C;
HinfSF.S = [HinfSF.C'*HinfSF.D1 HinfSF.C'*HinfSF.D2];

%% Gamma Search Algorithm

gamma_max = 20;
gamma_min = 1;
gamma = gamma_max;
to_test = 1;

while(abs(gamma_max - gamma_min)>1e-10)
    if to_test == 1
        gamma_max = gamma;
    else
        gamma_min = gamma;
    end
    gamma = (gamma_max + gamma_min)/2;   
    disp(['gamma_max = ',num2str(gamma_max,  '%12.10g'),'  gamma = ',num2str(gamma,  '%12.10g'),...
          '  gamma_min = ',num2str(gamma_min,  '%12.10g')])
    
    % Set Up Matrices
    HinfSF.R = [HinfSF.D1'*HinfSF.D1 HinfSF.D1'*HinfSF.D2;...
                HinfSF.D2'*HinfSF.D1 HinfSF.D2'*HinfSF.D2 - gamma^2];
    
    % Solve Riccati
    [P,CL_eig,K] = care(HinfSF.A,HinfSF.B_hat,HinfSF.Q,HinfSF.R,HinfSF.S);
    
    pos_test = 1;
    ev = eig(P);
    for i = 1:6
        if ev(i) < 0
            pos_test = 0;
        end
    end
    
    % Gain Matrix
    HinfSF.Kxopt = -[1 0]*HinfSF.R^-1*(HinfSF.B_hat'*P + HinfSF.S');
    
    HinfSF.Aref = HinfSF.A + HinfSF.B*HinfSF.Kxopt;
    
    % Check Valid Design
    eig_test = 1;
    ev_cl = eig(HinfSF.Aref);
    for i = 1:6
        if ev_cl(i) >= 0
            eig_test = 0;
        end
    end
    to_test = pos_test*eig_test;    
end

HinfSF.Kxref = HinfSF.Kxopt;
disp(['Gamma optimal = ',num2str(gamma,  '%12.10g')])
disp(['K optimal = ',num2str(HinfSF.Kxopt)])
Niter = 0;

EE=HinfSF.A'*P+P*HinfSF.A - (P*HinfSF.B_hat+HinfSF.S)*inv(HinfSF.R)*(HinfSF.B_hat'*P+HinfSF.S')+HinfSF.Q;
disp(['Norm EE opt = ',num2str(norm(EE))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(HinfSF.Kxref(2)>1 & Niter < 100)
    disp(['gamma = ',num2str(gamma),'  Kq = ',num2str(HinfSF.Kxref(2)),...
          '  # of iterations = ', num2str(Niter)])
    
    HinfSF.R = [HinfSF.D1'*HinfSF.D1 HinfSF.D1'*HinfSF.D2;...
                HinfSF.D2'*HinfSF.D1 HinfSF.D2'*HinfSF.D2 - gamma^2];
    
    [P,CL_eig,K] = icare(HinfSF.A,HinfSF.B_hat,HinfSF.Q,HinfSF.R,HinfSF.S);
    
    HinfSF.Kxref = -[1 0]*HinfSF.R^-1*(HinfSF.B_hat'*P + HinfSF.S');
    HinfSF.Aref = HinfSF.A + HinfSF.B*HinfSF.Kxref;

    gamma = gamma + 0.00005;
    Niter = Niter + 1;
end

disp(['Gamma = ',num2str(gamma, '%12.10g')])
disp(['K ref = ',num2str(HinfSF.Kxref)])

%% Closed-Loop System

HinfSF.Bref = HinfSF.E;
HinfSF.Cref = eye(6,6);
HinfSF.Dref = zeros(6,1);
HinfSF.sys_cl = ss(HinfSF.Aref,HinfSF.Bref,HinfSF.Cref,HinfSF.Dref);
disp(['Closed-Loop System Eigenvalues'])
damp(eig(HinfSF.Aref))

% Compute Riccati Eq solution

HinfSF.A,HinfSF.B_hat,HinfSF.Q,HinfSF.R,HinfSF.S

%A'X + 'XA - ('XB + S)R  (B'X + S') + Q = 0
EE=HinfSF.A'*P+P*HinfSF.A - (P*HinfSF.B_hat+HinfSF.S)*inv(HinfSF.R)*(HinfSF.B_hat'*P+HinfSF.S')+HinfSF.Q;
disp(['Norm EE = ',num2str(norm(EE))]);

%% Plotting Commands

t_step = [0:0.01:10]';
u_step = zeros(size(t_step));
u_step(201:300) = 32.174*ones(100,1);
u_step(601:800) = 3*32.174*ones(200,1);
[y,t] = lsim(HinfSF.sys_cl,u_step,t_step);
figure
plot(t,u_step,t,y(:,1),'LineWidth',2)
ylabel('Azc vs Az ~ fps2')
xlabel('Time ~ sec')
legend('Azc','Az')
grid on
if(save_plots == 1) saveppt2(plot_file_name); end   
figure
plot(t,y(:,2)*180/pi,'LineWidth',2)
ylabel('q ~ dps')
xlabel('Time ~ sec')
grid on
if(save_plots == 1) saveppt2(plot_file_name); end   
figure
plot(t,y(:,3),'LineWidth',2)
ylabel('dele ~ deg')
xlabel('Time ~ sec')
grid on
if(save_plots == 1) saveppt2(plot_file_name); end   
figure
plot(t,y(:,4),'LineWidth',2)
ylabel('dele dot ~ dps')
xlabel('Time ~ sec')
grid on
if(save_plots == 1) saveppt2(plot_file_name); end   
figure
plot(t,y(:,5),'LineWidth',2)
ylabel('xs')
xlabel('Time ~ sec')
grid on
if(save_plots == 1) saveppt2(plot_file_name); end   
figure
plot(t,y(:,6),'LineWidth',2)
ylabel('xt')
xlabel('Time ~ sec')
grid on
if(save_plots == 1) saveppt2(plot_file_name); end   

%*******************************************************
% General Form 
%*******************************************************  
%Close the loop to test the model
% Plant form  xdot = Apx + Bpu; 
%                      y = Cpx +Dpu
% Controller xcdot = Acxc + Bc1y + Bc2r
%                      u = Ccxc + Dc1y + Dc2r

  Ac = [ Ws.A       0.;
               0.    Wt.A]
  Bc1 = [ Ws.B*[1 0 0 0]
              Wt.B*[1 0 0 0]]
  Bc2 = [-Ws.B
                 0.]
  Cc   = [ HinfSF.Kxref(:,5:6) ]
  Dc1 = [ HinfSF.Kxref(:,1:4)]
  Dc2 = [ 0 ]
  
  
  Ap=Plant.Ap;
  Bp=Plant.Bp;
  Cp=eye(4);
  Dp=0*Cp*Bp;
  
%   w=logspace(-1,4,500);
  rtd=180/pi;

  
Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
Acl = [     (Ap+Bp*Z*Dc1*Cp)              (Bp*Z*Cc);
          (Bc1*(Cp+Dp*Z*Dc1*Cp))  (Ac+Bc1*Dp*Z*Cc)];
Bcl = [       Bp*Z*Dc2;
           (Bc2+Bc1*Dp*Z*Dc2)];
Ccl=[eye(4) zeros(4,2)];

Acl_sof_n = Acl;

w = logspace(-3,4,500);

t=linspace(0,2.0);
% Time Domain Analysis
[yy,x] = step(Acl,Bcl,Ccl,0.*Ccl*Bcl,1,t);
  t_hinf = t;
  y_hinf = yy(:,1);
  az = y_hinf;
  aze = abs(ones(size(az))-az);
  tr = 0.; ts= 0.;
  fv = aze(numel(aze));
  e_n = aze - fv*ones(size(aze)) - 0.36*ones(size(aze));
  tr = crosst(e_n,t);
  e_n = aze - fv*ones(size(aze)) - 0.05*ones(size(aze));
  ts = crosst(e_n,t);
  azmin = abs(min(az));
  azmax = abs(max(az));
  dmax = rtd*max(abs(x(:,3)));
  ddmax = rtd*max(abs(x(:,4)));
  data = [tr ts azmin azmax dmax ddmax];
% Frequency Domain Analysis
for ii=1:numel(w),
    GG = Cp*inv(sqrt(-1)*w(ii)*eye(size(Ap))-Ap)*Bp+Dp;
    KK = -Cc*inv(sqrt(-1)*w(ii)*eye(size(Ac))-Ac)*Bc1-Dc1;
    l(ii) = KK*GG;
  end;
  re_hinf = real(l);
  im_hinf = imag(l);
  rd = ones(size(l))+l;
  sr = ones(size(l))+ones(size(l))./l;
  rd_hinf = 20.*log10(abs(rd));
  sr_hinf = 20.*log10(abs(sr));
  rdmin = min(abs(rd));
  srmin = min(abs(sr));
  madb_hinf = 20.*log10(abs(l))';  
  pha_hinf = rtd*phase(l);
  wc = crosst(madb_hinf,w);
  data_hinf = [ data rdmin srmin  wc];
  

figure
plot(t,yy(:,1),'b','LineWidth',2);grid
xlabel('Time (sec)');
ylabel('Az (fps2)');
if(save_plots == 1) saveppt2(plot_file_name); end   

figure
plot(t,yy(:,2)*rtd,'b','LineWidth',2);grid
xlabel('Time (sec)');
ylabel('Pitch Rate (dps)');
if(save_plots == 1) saveppt2(plot_file_name); end   

figure
plot(t,yy(:,3)*rtd,'b','LineWidth',2);grid
xlabel('Time (sec)');
ylabel('Elevon (deg)');
if(save_plots == 1) saveppt2(plot_file_name); end   

figure
plot(t,yy(:,4)*rtd,'b','LineWidth',2);grid
xlabel('Time (sec)');
ylabel('Elevon Rate(dps)');
if(save_plots == 1) saveppt2(plot_file_name); end   


