function Y = rbf2mat_v2(X,T,W)
% -----------------------------------------------------------
% Purpose: RBF2MAT_V2 computes RBF functions.
%
% Call:  Y = rbf2mat_v2(X,T,W)
% 
% Inputs:
%   X - (nxN) - matrix of N n-dim patterns.
%   T - (nxm) - matrix of m n-dim RBF centers.
%   W - (nxn) - positive definite weight matrix.
%
% Outputs:
%   Y - (Nxm) - matrix of N rows of m rbf-function values.
%
% Created : Eugene Lavretsky, 03-12-02
%
% Modified: 01-16-03, Eugene Lavretsky: RBF generic definition
% -----------------------------------------------------------

[nX,mX] = size(X);
[nT,mT] = size(T);

if nX~=nT,
  error('# of rows in X and T must be the same');
else   
  N = mX;
  m = mT;
end;

Y = zeros(N,m);	% output array pre-allocation

for i=1:N,
  for j = 1:m,
    z = X(:,i)-T(:,j);
    Y(i,j) = exp(-z'*W*z);
  end  
end;

return