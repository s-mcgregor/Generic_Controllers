close all;
clear

%Robust Control - HW2, Problem 1
%*************************************************************************
% Plant Model
%*************************************************************************
% State Names
% -----------
% AZ fps2
% q rps
% Dele deg
% Dele dot dps
% Input Names
% -----------
% Dele cmd deg
Ap = [ -0.576007 -3255.07 4.88557 9.25796;
    -0.0410072 -0.488642 -2.03681 0 ;
    0 0 0 1 ;
    0 0 -8882.64 -133.266 ];
Bp = [ 0 ; 0 ; 0 ; 8882.64];
Cp = [ 1 0 0 0; 
    0 1 0 0];
Dp = 0.*Cp*Bp;

[nCp,nAp] = size(Cp);
%*******************************************************
% Static Output Feedback With New Form To Test Prior To Substituting
%The Big Plant Model
%*******************************************************
%Close the loop to test the model
% Plant form xdot = Apx + Bpu; y = Cpx +Dpu
% Controller xcdot = Acxc + Bc1y + Bc2r
% u = Ccxc + Dc1y + Dc2r
Ac = [ 0 ];
Bc1 = [ -1 0];
Bc2 = [ 1 ];
Cc = [ 0.0107349];
Dc1 = [ -0.0411729 11.4003];
Dc2 = [ 0 ];

% SS model of loop gain L at the plant output
zero_vector_Ly = zeros(size(Ac,1),size(Ap,1)); %Create zero vector
Aout = [Ap Bp*Cc;  zero_vector_Ly Ac];
Bout = [Bp*Dc1; Bc1];
Cout = -[Cp Dp*Cc]; %-
Dout = -[Dp*Dc1]; %-

% SS model of loop gain Lu at the plant input
zero_vector_Lu = zeros(size(Ap,1),size(Ac,1)); %Create zero vector
Ain = [ Ap zero_vector_Lu;   Bc1*Cp Ac];
Bin = [ Bp; Bc1*Dp];
Cin = -[ Dc1*Cp Cc]; %-
Din = -[Dc1*Dp]; %-

% Lu & Ly transfer functions, sensitivity, and complimentary sensitivity
L_input = ss(Ain,Bin,Cin,Din); % Lu
S_input = inv(eye(size(L_input))+L_input); % Sensitivity
T_input = eye(size(L_input)) - S_input; % Complimentary Sensitivity

L_output = ss(Aout,Bout,Cout,Dout); % Ly
S_output = inv(eye(size(L_output))+L_output); % Sensitivity
T_output = eye(size(L_output)) - S_output; % Complimentary Sensitivity

% 1a Break loop at plant input (Lu); Plot Bode Plot, id gain and phase margins, LGXF, PXF
margin(L_input)
[GM, PM_deg, wc_GM, wc_Pm] = margin(L_input);

% 1b Plot Nyquist of Lu, id gain and phase margin, loop gain and phase xover
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
RDu_min = min(min(abs(RDu_HS)))
SRu_min = min(min(abs(SRu_HS)))

RDu_nGM = 1/(1+RDu_min);
RDu_pGM = 1/(1-RDu_min);
RDu_Pha = 2*asin(RDu_min/2);
RDu_nGM_dB = 20*log10(RDu_nGM)
RDu_pGM_dB = 20*log10(RDu_pGM)
RDu_Pha_deg = 180*RDu_Pha/pi 

SRu_nGM = 1-SRu_min;
SRu_pGM = 1+SRu_min;
SRu_Pha = 2*asin(SRu_min/2);
SRu_nGM_dB = 20*log10(SRu_nGM)
SRu_pGM_dB = 20*log10(SRu_pGM)
SRu_Pha_deg = 180*SRu_Pha/pi 
disp('  ')

%disp('Classical Margins')
%allmargin(sys_u)


%Sigma(I+Lu)
figure
sigma(eye(size(L_input))+L_input);
title('Sigma (I + Lu)')
legend([' min(I+Lu) = ' num2str(RDu_min)],'Location','Best');

%Sigma(I+Inv(Lu))
figure
sigma(eye(size(L_input))+inv(L_input));
title('Sigma (I + inv(Lu))')
legend([' min(I+invLu) = ' num2str(SRu_min)],'Location','Best');
