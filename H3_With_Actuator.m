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

%SS model of loop gain at the plant input
Ain = [ Ap 0.*Bp*Cc;  Bc1*Cp Ac];
Bin = [ Bp; Bc1*Dp];
Cin = -[ Dc1*Cp Cc];%change sign for loop gain
Din = -[ Dc1*Dp];
Sys_Lu = ss(Ain, Bin, Cin, Din);

[reLu,imLu]   = nyquist(Ain,Bin,Cin,Din,1,w);
[magLu,phaLu] = bode(Ain,Bin,Cin,Din,1,w);


%SS model of loop gain at the plant output
Aout = [ Ap Bp*Cc;  0.*Bc1*Cp Ac];
Bout = [ Bp*Dc1; Bc1];
Cout = -[ Cp Dp*Cc];%change sign for loop gain
Dout = -[ Dp*Dc1];
Sys_Ly = ss(Aout, Bout, Cout, Dout);

% Plant  Freqeuncy Domain Analysis using looping method which evaluates the
% frequency response at each frequency s = j*w
for i=1:numel(w),
    s = sqrt(-1)*w(i);
    GG = Cp*inv(s*eye(size(Ap))-Ap)*Bp+Dp; % Plant
    KK = Cc*inv(s*eye(size(Ac))-Ac)*Bc1+Dc1; % Controller
    
    Lu(i) = -KK*GG; % This should produce the same result as Sys_Lu above
    RDu(i) = 1.+Lu(i);
    SRu(i) = 1.+1./Lu(i);
    
    Ly = Cout*inv(s*eye(size(Aout))-Aout)*Bout+Dout;
    Sy = inv(eye(size(Ly))+Ly);
    Ty = Ly*Sy;
    
    Sy_1(i) = abs(Sy(1,1));
    Ty_1(i) = abs(Ty(1,1));
    
    Ly2 = -GG*KK; % This should produce the same result as Sys_Ly, Ly above
    Rdy2 = eye(size(Ly2))+Ly2;
    RDy2_min(i) = min(svd(Rdy2));
    Sy2 = inv(Rdy2);
    Ty2 = Ly2*Sy2;
    Sy_2(i) = abs(Sy2(1,1));
    Ty_2(i) = abs(Ty2(1,1));
end

 T  = freqresp(sys_cl,w); % Complementary Sensitivity
 S = 1 - T; % Sensitivity



magdB = 20.*log10(magLu);
wc = crosst(magdB,w);
rdu_min = min(abs(RDu));
sru_min = min(abs(SRu));
RDy_out_min = min(RDy2_min);
S_out_max = max(Sy_1);
T_out_max = max(Ty_1);
S2_out_max = max(Sy_2);
T2_out_max = max(Ty_2);

%Compute closed loop eigenvalues and vectors
%disp('Closed loop eigenvalues and vectors')
%[V,D] = eig(Acl);
%disp(['Closed Loop Eigenvalues =    ' num2str(D)])
%disp(['Closed Loop Eigenvectors = ' num2str(V)])

figure('Name','Bode Using Margin Command'),
%Compute classical margins from the bode plot
disp('Classical margins using Matlab command')
margin(Sys_Lu)

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
    num2str(pm)  ' deg' ])


% Plot
%
Tit1 = 'Problem 1 Frequency Response';
ngm = -1/neg_gm;
pgm = -1/pos_gm;

figure('Name','Nyquist Plot'),
plot(xx1,yy1,'k:',real(Lu),imag(Lu),'k',reLu,imLu,'b',...
    [ngm -1.],[0. 0.],'r',[-1 pgm],[0. 0.],'c','LineWidth',2);grid
axis([-3 3 -3 3]);
xlabel('Re(Lu)')
ylabel('Im(Lu)')
legend('Unit Circle at -1,j0','Looping method','SS method','Neg SV Margin','Pos SV Margin');
title(Tit1)

figure('Name','Bode Magnitude'),
semilogx(w,magdB,'b','LineWidth',2);grid
legend(['wc = ' num2str(wc) ' rps'],'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag (dB)')
title(Tit1)


figure('Name','Bode Phase'),
semilogx(w,phaLu,'b','LineWidth',2);grid
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
title(Tit1)


figure('Name','I+Lu'),
semilogx(w,20*log10(abs(RDu)),'b','LineWidth',2);grid
legend(['min I+Lu = ' num2str(rdu_min)],'Location','Best');
xlabel('Frequency (rps)')
ylabel('RDu Mag dB')
title(Tit1)

figure('Name','I+invLu'),
semilogx(w,20*log10(abs(SRu)),'b','LineWidth',2);grid
legend(['min I+invLu = ' num2str(sru_min)],'Location','Best');
xlabel('Frequency (rps)')
ylabel('SRu Mag dB')

figure('Name','I+Ly'),
semilogx(w,20*log10(abs(RDy2_min)),'b','LineWidth',2);grid
legend(['min I+Ly = ' num2str(RDy_out_min)],'Location','Best');
xlabel('Frequency (rps)')
ylabel('RDy Mag dB')

figure('Name','Comp Sensitivity T=Y/R')
T_az_max = max(abs(squeeze(T(1,1,:))));
semilogx(w,20*log10(abs(squeeze(T(1,1,:)))),'LineWidth',2);grid
tex1 = ['T-Az inf norm = ' num2str(T_az_max)];
legend(tex1,'Location','Best');
xlabel('Frequency (rps)')
ylabel('Co-Sens Mag (dB)')
title(Tit1)

figure('Name',' Sensitivity S=E/R')
S_az_max = max(abs(squeeze(S(1,1,:))));
semilogx(w,20*log10(abs(squeeze(S(1,1,:)))),'LineWidth',2);grid
tex1 = ['S-Az inf norm = ' num2str(S_az_max)];
legend(tex1,'Location','Best');
xlabel('Frequency (rps)')
ylabel('Sens Mag (dB)')
title(Tit1)
