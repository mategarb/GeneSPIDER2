function [hatTheta hatXi hatPhi] = tls(Xi,Phi)
% hatTheta = tls(Xi,Phi)
% Basic total least squares estimation algorithm (cite:Markovsky2007)
% 
% From the model $(\Phi - \Upsilon) \times \Theta = (\Xi - \Pi)$
% where
% $\Xi = -P^T$
% $\Phi = Y^T$
% $\Theta = A^T$
%

m        = size(Phi,1);        % number of x,y data pairs
X        = Xi;
Y        = Phi;
n        = size(Xi,2);         % n is the width of X (X is m by n)
Z        = [Xi Phi];           % Z is X augmented with Y.
[U S V]  = svd(Z,0);           % find the SVD of Z.
VXY      = V(1:n,1+n:end);     % Take the block of V consisting of the first n rows and the n+1 to last column
VYY      = V(1+n:end,1+n:end); % Take the bottom-right block of V.
hatTheta = -VXY/VYY;

S1 = S;
S2 = S;

S1(n+1:end,n+1:end) = 0;
% S2(1:n,1:n) = 0;

hatC = U*S1*V;

hatXi = hatC(:,1:n);
hatPhi = hatC(:,n+1:end);
