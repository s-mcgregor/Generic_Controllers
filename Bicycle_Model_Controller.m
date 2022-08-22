%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dynamic Model of Ackermann Steering %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% Define required vehicle parameters 

mass = 100; % units
Ca_F = 1; % cornering stiffness, front ; units
Ca_R = 1; % cornering stiffness, rear ; units
L_F = 10; % CG to front axle distance ; units
L_R = 10; % CG to rear axle distance ; units
I_z = 1000; % vehicle Moment of Interia (z) ; units

% Vehicle speed ot be modified from two extremes
Vx = 1; % units

% Define state equations
% State 1: e1: distance from the CG to center line of lane
% State 2: e2: orientation error of vehicle w.r.t road
% State 3: e1_dot
% State 4: e2_dot
% Input 1: Front Wheel Angle (d)
% Input 2: Desired Vehicle Orientation (psi)

Ap = [0 1 0 0;
    0, -((2*Ca_F)+(2*Ca_R))/(mass*Vx), ((2*Ca_F)+(2*Ca_R))/mass, ((-2*Ca_F*L_F)+(2*Ca_R*L_R))/(mass*Vx);
    0 0 0 1;
    0, -((2*Ca_F*L_F)-(2*Ca_R*L_R))/(I_z*Vx), ((2*Ca_F*L_F)-(2*Ca_R*L_R))/I_z, -((2*Ca_F*L_F*L_F)+(2*Ca_R*L_R*L_R))/(I_z*Vx)];
B_d = [0; (2*Ca_F)/mass ; 0 ; (2*Ca_F*L_F)/I_z];
B_psi = [0; (-((2*Ca_F*L_F)-(2*Ca_R*L_R))/(mass*Vx)) - Vx; 0; -((2*Ca_F*L_F*L_F)+(2*Ca_R*L_R*L_R))/(I_z*Vx)];
Bp = [B_d B_psi];
[nx,~] = size(Ap);
Cp = eye(nx); % start with full state output
Dp = 0.*Cp*Bp;

%Check Controlability
if nx == rank(ctrb(Ap,Bp))
    disp('Controllable')
else
    disp('Not Controllable')
end

%Check Observability
if nx == rank(obsv(Ap,Cp))
    disp('Observable')
else
    disp('Not Observable')
end

% Setup range of penalties for the LQR
Q=0.*Ap; %sets up a matrix Q that is the size of Aw and is all 0s
Q(1,1) = 1;
Q(2,2) = 1;
R=eye(2);
R(1,1) = 100;
[Kx_lqr,~,~]=lqr(Ap,Bp,Q,R);

 % populate the controller matrices
    Ac = zeros(2);
    Bc1 = [0. 0. 0. 0.;
           0. 0. 0. 0.];
    Bc2 = 0.;
    Cc = zeros(2);
    Dc1 = [-Kx_lqr];
    Dc2 = [1];

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
    
    L_input = ss(A_Lu,B_Lu,C_Lu,D_Lu); % Loop Gain at Plant Input
    S_input = inv(eye(size(L_input))+L_input); % Sensitivity
    T_input = eye(size(L_input)) - S_input; % Complimentary Sensitivity

% Stability margin information    
% [GM, PM_deg, wc_GM, wc_Pm] = margin(L_input)
% Margin Requires a SISO system

% Plot Nyquist of Lu
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