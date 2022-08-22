clc
clear all
close all

% ME598 Robust Control - Homework 3, Problem 1

Za_V = -1.3046;
Ma   = 47.711;
Zd_V = -.2142;
Md   = -104.83;
V    = 886.78;
Za   = V*Za_V;
Zd   = V*Zd_V;
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

% Form Squiggle Matrix

zeros_A1 = zeros(size(Cp,1));
zeros_A2 = transpose(zeros(size(Cp)));

% z = [e; Az_dot; q_dot; de_dot de_ddot]'

%A_sq = [zeros_A1 Cp ; zeros_A2 Ap];
%B_sq = [Dp; Bp];

zero_vector = [0 0 0 0]';
A_sq = [0 1 0 0 0; zero_vector Ap]
B_sq = [0;Bp];
C_sq = eye(size(A_sq)); %output all states
D_sq = 0.*C_sq*B_sq;

F = [-1 ; zeros(size(A_sq,1)-1,1)]; % [-1; 0; ...; 0]

% Objective is to create design to track acceleration. Start with R = 1 and
% selecting a Q matrix that penalizes the error state e
% Performance index J = integral (transpose(z)*Q*z + mu^2)dTao
% Start by setting Q11 = something, and all other elements = 0. Sweep q11
% from large to small and examine the closed loop properties.

Q = zeros(size(A_sq));
q11 = 1;
Q(1,1) = q11;
R= [1];

%Check Observability
nx = size(A_sq,1);
if nx == rank(obsv(A_sq,Q^.5))
    disp('Observable')
else
    disp('Not Observable')
end

sysSq = ss(A_sq,B_sq,C_sq,D_sq);
K = lqr(sysSq,Q,R);
sysCL = ss(A_sq-B_sq*K,[B_sq],C_sq,D_sq);


