function varargout = cvxcls(varargin)
% mle does cvx optimization on a network matrix
%
%   Input Arguments: mle(data,initA)
%   ================
%   data   = The data in a struct array. Must contain the
%            covariance of 'nY' and 'P' as 'cvY' and 'cvP'
%   initA  = Precomputed network matrix as cell array.
%
%   Output Arguments: Aest
%   =================
%   0 or 1: if only data and required are given then an if
%   statement can be issued on the scripts to see if the method can
%   be used with the data.
%   Aest = The estimated network matrices as cell array.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse input arguments %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawZeta = 0;
net = [];
for i=1:nargin
    if isa(varargin{i},'GeneSpider.Dataset')
        data = varargin{i};
    elseif isa(varargin{i},'GeneSpider.Network')
        net = varargin{i};
    else
        initA = varargin{i};
    end
end

if ~exist('data')
    error('needs a data set')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parse data to use
%% Parse data to use
P = data.P;
Y = response(data,net);

[cvY,cvP] = cov(data);

if ~isempty(data.cvY)
    cvY = data.cvY;
end

if ~isempty(data.cvP)
    cvP = data.cvP;
end

%% run

tol = eps;
[nGenes, nExp] = size(data.P);

if ~exist('initA','var') && ~isempty(net)
    initA = lsqr(data,net);
end
e = 1e-6;

for i=1:size(initA,3)
    cvx_clear
    try
        Ahat = initA(:,:,i);
        [H, D, H] = svd(Ahat * cvY * Ahat' + cvP + e*eye(size(Ahat)));
        dD = diag(D);
        Correction = H * sqrt(diag( 1 ./dD ));
        S = abs(Ahat) < tol; % Constraints elements to be 0 to avoid bias term.
        cvx_begin quiet
        cvx_precision(eps*100)
        variables Ahat(nGenes, nGenes) eta(nGenes, nExp)
        minimize norm(eta'*Correction,'fro') % + 1e-12.*sum(sum(abs(eta)))
        subject to
        Ahat(S) == 0 % Forcing elements to zero
        Ahat*Y+P == eta
        cvx_end
    catch ME
        rethrow(ME)
        % warning('%s\nPure MLE failed. Introducing weak lasso constraint.',ME.message)
        Ahat = initA(:,:,i);
        [H, D, H] = svd(Ahat * cvY * Ahat' + cvP + e*eye(size(Ahat)));
        dD = diag(D);
        Correction = H * sqrt(diag( 1 ./dD ));
        S = abs(Ahat) < tol; % Constraints elements to be 0 to avoid bias term.
        fprintf('%s\n',ME.message)
        cvx_clear % Clears the CVX environment
        cvx_begin quiet
        variables Ahat(nGenes, nGenes) eta(nGenes, nExp)
        minimize norm(eta'*Correction,'fro') + 1e-15.*sum(sum(abs(eta)))
        subject to
        Ahat(S) == 0 % Forcing elements to zero
        Ahat*Y+P == eta
        cvx_end
    end
    Aest(:,:,i) = full(double(Ahat));
end

cvx_clear
varargout{1} = Aest;