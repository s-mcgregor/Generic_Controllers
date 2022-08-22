function [KX_hat_dot, Theta_hat_dot, ...
          KX_hat_new, Theta_hat_new, u_ad] = ...
          ofmrac(Y, Y_hat, x_hat, u_bl, ...
                           KX_hat, Theta_hat,  ...
                           KX_hat_dot_p, Theta_hat_dot_p,  ...
                           RBF_centers, RBF_W, MixY, ...
                           Gamma_X, Gamma_Theta, ...
                           gamma_X, gamma_Theta, dead_zone, ...
                           KX_max, Theta_max, dt)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Purpose : Output Feedback MRAC with projection, e-modification, and dead-zone
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Function call:
%       [KX_hat_dot, Theta_hat_dot, ...
%        KX_hat_new, Theta_hat_new, u_ad] = ...
%           ofmrac(Y, Y_hat, X_hat, X_ref, ...
%                            KX_hat, Theta_hat,  ...
%                            KX_hat_dot_p, Theta_hat_dot_p,  ...
%                            RBF_centers, RBF_W, MixY, ...
%                            Gamma_X, Gamma_Theta, ...
%                            gamma_X, gamma_Theta, dead_zone, ...
%                            KX_max, Theta_max, dt)
%
% Inputs  :
%   Y                   - (nYx1)  output
%   Y_hat               - (nYx1)  estimated output
%   x_hat               - (nXx1)  estimated state
%   u_bl                - (mx1)   baseline control input
%   KX_hat              - (nXxm)  X feedback gain matrix,
%   Theta_hat           - (Nxm)   matrix of estimated parameters for K0 uncertainty
%   KX_hat_dot_p        - (nXxm)  matrix of previous KX_hat time derivative
%   Theta_hat_dot_p     - (Nxm)   matrix of previous Theta_hat time derivative
%   RBF_centers         - (nXxN)  matrix of N n-dim RBF centers
%   RBF_W               - (NxN)   RBF positive definite weight matrix of sigmas
%   MixY                - (nYxnY) output mixing matrix for adaptive laws
%   Gamma_X             - (nXxnX) positive definite matrix of learning rates for KX_hat adaptation
%   Gamma_Theta         - (NxN)   positive definite matrix of learning rates for Theta_hat adaptation
%   gamma_X             - (1x1)   e-mod rate for KX_hat adaptation
%   gamma_Theta         - (1x1)   e-mod rate for Theta_hat adaptation
%   dead_zone           - (1x1)   dead-zone tolerance
%   KX_max              - (1x1)   KX_hat norm upper bound
%   Theta_max           - (1x1)   Theta_hat norm upper bound
%   dt                  - (1x1)   simulation step time, sec
%
% Outputs :
%   KX_hat_dot          - (nXxm)  matrix of current KX_hat time derivative
%   Theta_hat_dot       - (Nxm)   matrix of current Theta_hat time derivative
%   KX_hat_new          - (nXxm)  next X feedback gain matrix
%   Theta_hat_new       - (Nxm)   next Theta_hat matrix
%   u_ad                - (mx1)   vector of adaptive control inputs
%
% Created : 11-02-09, Eugene Lavretsky
%
% Modified:
%   10-22-11, Eugene Lavretsky
%   1) Changed e to ebar
%
%   11-06-11, Eugene Lavretsky
%   1) added return to if-else logic
%   2) changed X_ref to u_bl
%
%   01-22-12, Eugene Lavretsky
%   Changed Projection Operator (moved Gamma inside)
%   Theta_dot = Proj(Theta, Gammma*Theta);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% RBF regressor vector
RBF = rbf2mat_v2(x_hat, RBF_centers, RBF_W);

% add bias
PSI = [RBF, 1]';

% RBF approximation
F_hat = Theta_hat'*PSI;

% Linear regressor vector
X_hat = [u_bl; x_hat];

% Adaptive control
u_ad = -KX_hat'*X_hat - F_hat;

m = length(u_ad);  % number of controls

% Output estimation error
eY = Y - Y_hat;

% Training error and its magnutude
ebar      = eY'*MixY;
norm_ebar = norm(ebar);

% Projection epsilon, (defines smooth transition region)
proj_eps = 0.1; %  10 percent

% Output feedback adaptive laws
if norm_ebar <= dead_zone,
  KX_hat_dot    = KX_hat_dot_p*0;
  Theta_hat_dot = Theta_hat_dot_p*0;

  KX_hat_new    = KX_hat;
  Theta_hat_new = Theta_hat;
  return
else
  KX_hat_dot    =  Gamma_X*(X_hat*ebar - gamma_X*norm_ebar*KX_hat);
  Theta_hat_dot =  Gamma_Theta*(PSI  *ebar - gamma_Theta*norm_ebar*Theta_hat);
  
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % Parameter Projection logic starts here
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for j = 1:m,
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % KX_hat Projection
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    theta    = KX_hat(:,j);
    theta_max= KX_max(j);
    proj_tol = proj_eps*theta_max^2;                         % projection tolerance
    
    y        = KX_hat_dot(:,j);
    
    f        = ((1+proj_eps)*(theta')*(theta) - theta_max^2)/proj_tol;
    grad_f   = 2*(1+proj_eps)*theta/proj_tol;
    
    if f > 0 && y'*grad_f > 0,
      proj_y = y - Gamma_X*(grad_f)*(grad_f')*y*f/(grad_f'*Gamma_X*grad_f);
      KX_hat_dot(:,j) = proj_y;
    end;
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Theta_hat Projection
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    theta    = Theta_hat(:,j);
    theta_max= Theta_max(j);
    proj_tol = proj_eps*theta_max^2;                         % projection tolerance
    
    y        = Theta_hat_dot(:,j);
    
    f        = ((1+proj_eps)*(theta')*(theta) - theta_max^2)/proj_tol;
    grad_f   = 2*(1+proj_eps)*theta/proj_tol;
    
    if f > 0 && y'*grad_f > 0,
      proj_y = y - Gamma_Theta*(grad_f)*(grad_f')*y*f/(grad_f'*Gamma_Theta*grad_f);
      Theta_hat_dot(:,j) = proj_y;
    end;
    
  end;
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % Projection logic ends here
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end;

% AB-2 one-step-forward integration
KX_hat_new    = KX_hat    + dt*(1.5*KX_hat_dot - 0.5*KX_hat_dot_p);
Theta_hat_new = Theta_hat + dt*(1.5*Theta_hat_dot - 0.5*Theta_hat_dot_p);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% protection against numerical integration errors
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sqrt_1p1 = sqrt(1 + proj_eps);
for j = 1:m,
  norm_KX_j = norm(KX_hat_new(:,j));
  if norm_KX_j > sqrt_1p1*KX_max(j),
    KX_hat_new(:,j) = KX_hat_new(:,j) / norm_KX_j * sqrt_1p1*KX_max(j);
  end;
  
  norm_Theta_j = norm(Theta_hat_new(:,j));
  if norm_Theta_j > sqrt_1p1*Theta_max(j),
    Theta_hat_new(:,j) = Theta_hat_new(:,j) / norm_Theta_j * sqrt_1p1*Theta_max(j);
  end;
end;

return;