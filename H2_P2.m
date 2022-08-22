clear all
close all

num = [-1];
den = [1 1];
L_input = tf(num,den); % Lu
S_input = inv(eye(size(L_input))+L_input); % Sensitivity
T_input = eye(size(L_input)) - S_input; % Complimentary Sensitivity
Neg_T_input =-1*T_input;
sigma(inv(Neg_T_input))
