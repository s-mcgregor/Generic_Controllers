num1 = [0.5 1];
den1 = [1 2.5 1];
num2 = [0.5 1];
den2 = [1 2.5 1 0];

tf(num1,den1)
[A,B,C,D] = tf2ss(num1,den1);

%[num,den] = ss2tf(Acl,Bcl,Ccl,Dcl);
%H = tf(num,den);
%bode(H)
%margin(H)
%nyquist(H)

L = ss(A,B,C,D); % L
S = inv(eye(size(L))+L); % Sensitivity
T = eye(size(L)) - S; % Complimentary Sensitivity
bodemag(L)
hold on
bodemag(S)
bodemag(T)
legend('|L|','|S|','|T|')
title('Bode Magnitude Diagram - Type 0')