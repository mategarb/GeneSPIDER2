function [hatTheta hatXi hatPhi] = tls(Phi,Xi)
% Basic total least squares estimation algorithm (Markovsky and Van Huffel 2007 http://www.sciencedirect.com/science/article/pii/S0165168407001405)
% hatTheta = tls(Xi,Phi)
% 
% From the model $$(\Phi - \Upsilon) \times \Theta = (\Xi - \Pi)$$
% where
% $\Xi = -P^T$
% $\Phi = Y^T$
% $\Theta = A^T$
%

m        = size(Xi,1);         % number of x,y data pairs
n        = size(Phi,2);        % n is the width of X (X is m by n)
Z        = [Phi Xi];           % Z is X augmented with Y.
[U S V]  = svd(Z,0);           % find the SVD of Z.
VXY      = V(1:n,1+n:end);     % Take the block of V consisting of the first n rows and the n+1 to last column
VYY      = V(1+n:end,1+n:end); % Take the bottom-right block of V.
hatTheta = -VXY/VYY;            % Generate TLS estimate

S1 = S;
S1(n+1:end,n+1:end) = 0;
hatC = U*S1*V;

hatXi = hatC(:,1:n);
hatPhi = hatC(:,n+1:end);


% S2 = S;
% S2(1:n,1:n) = 0;
% deltaC = U*S2*V;

