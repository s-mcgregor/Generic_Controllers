clear all
close all

% See Example 5.1 (Page 100/101 in Orange book)
% Given Values in homework problem
Ka = -0.0015;
Kq = -0.32;
az = 2.0;
aq = 6.0;
V = 886.78;
ZaV = -1.3046;
ZdV = -0.2142;
Ma = 47.7109;
Md = -104.8346;
% Calculated Values
Za = ZaV*V;
Zd = ZdV*V;

% Plant Model
% State 1: alpha- Angle of Attack (rad)
% State 2: q- Body Pitch Rate (rps)
% Output 1: Az- Acceleration (fps2)
% Output 2: q- Body Pitch Rate (rps)
Ap = [ZaV 1; Ma 0];
Bp = [ZdV; Md];
Cp = [Za 0; 0 1];
Dp = [Zd;0];

[nx,nx]=size(Ap); %number of states (nx)

% Common Controller Model
Ac = [0 0; Kq*aq 0];
Bc1 = -[Ka*az 0 ; Ka*Kq*aq Kq*aq];
Bc2 = [Ka*az ; Ka*Kq*aq];
Cc = [Kq 1];
Dc1 = -[Ka*Kq Kq];
Dc2 = [Ka*Kq];

% Close the loop
% Plant form  xdot = Apx + Bpu;
% y = Cpx +Dpu
% Controller xcdot = Acxc + Bc1y + Bc2r
% u = Ccxc + Dc1y + Dc2r
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
else
  disp('Closed-Loop System is stable');
  disp(['  Dominant Re(eig) = ', num2str(max(real(eig(Acl))))]);
end;

step(sys_cl)
[y,t]=step(sys_cl);  % t is time vector
y11=y(:,1);
figure
plot(t,y11);
title('Step Response - Output 1 - Acceleration')
xlabel('Time (s)') 
ylabel('Acceleration (fps2)') 

% SS model of loop gain Lu at the plant input
zero_vector_Lu = zeros(size(Ap,1),size(Ac,1)); %Create zero vector
Ain = [ Ap zero_vector_Lu;   Bc1*Cp Ac];
Bin = [ Bp; Bc1*Dp];
Cin = -[ Dc1*Cp Cc]; %-
Din = -[Dc1*Dp]; %-

% Lu transfer function, sensitivity, and complimentary sensitivity
L_input = ss(Ain,Bin,Cin,Din); % Lu
S_input = inv(eye(size(L_input))+L_input); % Sensitivity
T_input = eye(size(L_input)) - S_input; % Complimentary Sensitivity
Neg_T_input =-1*T_input;

% Lu Stability margins
figure
margin(L_input)
[GM, PM_deg, wc_GM, wc_Pm] = margin(L_input);

% Return Difference and SV margins
w = logspace(-1,3,500);

for i=1:numel(w),
    s = sqrt(-1)*w(i);
    GG = Cp*inv(s*eye(size(Ap))-Ap)*Bp+Dp;
    KK = Cc*inv(s*eye(size(Ac))-Ac)*Bc1+Dc1;
    Lu_HS(i)  = -KK*GG;
    RDu_HS(i)  = 1. + Lu_HS(i);
    SRu_HS(i) = 1. + 1./Lu_HS(i);
end

disp('  ')
disp('SV Margins')
RDu_min = min(min(abs(RDu_HS)));
SRu_min = min(min(abs(SRu_HS)));

RDu_nGM = 1/(1+RDu_min);
RDu_pGM = 1/(1-RDu_min);
RDu_Pha = 2*asin(RDu_min/2);
RDu_nGM_dB = 20*log10(RDu_nGM);
RDu_pGM_dB = 20*log10(RDu_pGM);
RDu_Pha_deg = 180*RDu_Pha/pi ;

SRu_nGM = 1-SRu_min;
SRu_pGM = 1+SRu_min;
SRu_Pha = 2*asin(SRu_min/2);
SRu_nGM_dB = 20*log10(SRu_nGM);
SRu_pGM_dB = 20*log10(SRu_pGM);
SRu_Pha_deg = 180*SRu_Pha/pi ;
disp('  ')


figure
sigma(inv(Neg_T_input))
hold on

tao = .01;
num = [-tao 0];
den = [tao 1];
D01 = tf(num,den);
sigma(D01);

tao = .05;
num = [-tao 0];
den = [tao 1];
D05 = tf(num,den);
sigma(D05);

tao = .1;
num = [-tao 0];
den = [tao 1];
D10 = tf(num,den);
sigma(D10);

legend({'1/sigma(M)','sigma(Delta); Tao = .01','sigma(Delta); Tao = .05','sigma(Delta); Tao = .10'},'Location','best')
