% MIMO RSLQR TEST
% Arbitrary Controllable/Observable 4 State Matrix
Ap = [20 30 40 50; 40 30 20 10; 0 -20 -50 10; 5 -5 5 -5];
Bp1 = [1; 10; 5; 10];
Bp2 = [5; 5; 5; 5];
Bp = [Bp1 Bp2];
Cp = eye(size(Ap));
Dp = 0.*Cp*Bp;

[nx,~] = size(Ap);
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

% Create wiggle system
% Lecture 10 - Optimal Control - Slide 37

C = [1 0 0 0; 0 1 0 0];
D = 0.*C*Bp;
zcp = zeros(size(C,1));
zcp(1,2)= 1; % I would have thought this to be zero since constant command
zap = zeros(size(Ap,1),2);
Aw = [zcp C;
      zap Ap];
Bw = [D; Bp];
F = [-1 ; zeros(size(Aw,1)-1,1)]; % [-1; 0; ...; 0]
Q = zeros(size(Aw,1));
Q(1,1) = 10;
Q(2,2) = 100;
R = 1;

[Kx_lqr,~,~]=lqr(Aw,Bw,Q,R);
Kx_lqr
K_X = Kx_lqr(:,3:6)
K_I = Kx_lqr(:,1:2)