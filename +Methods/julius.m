function varargout = julius(varargin)
% function estA = julius(data,require,zetavec,rawZeta)
% evalGLasso will evaluate the method glasso and produce an
% estimated network matrix.
%
%   Input Arguments: julius(data,net,zetavec [,rawZeta,'initA',initA])
%   ================
%   data:    GeneSpider.Dataset
%   net:     GeneSpider.Network
%   zetavec: method parameter for tuning the network fitting. (optinal)
%   rawZeta: logical to determine if the zeta values should be
%            converted.  default = false
%   initA  = Precomputed network matrix as 3d array. input as string value pair
%
%   Output Arguments: estA
%   =================
%   estA: the estimated networks as a 3d array.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse input arguments %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawZeta = 0;

for i=1:nargin
    if isa(varargin{i},'GeneSpider.Dataset')
        data = varargin{i};
    elseif isa(varargin{i},'GeneSpider.Network')
        net = varargin{i};
    elseif isa(varargin{i},'logical')
        rawZeta = varargin{i};
    elseif isa(varargin{i},'char')
        i = i+1;
        initA = varargin{i};
    else
        zetavec = varargin{i};
    end
end

if ~exist('data')
    error('needs a data set')
end
if ~exist('net')
    error('needs a network')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine how to handle zeta %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~rawZeta
    zetaRange = [];
    tol = 1e-6;
    zmax = 1;
    zmin = 0;
    % find zero network
    estA = Methods.julius(data,net,zmax,logical(1));
    while nnz(estA) > 0 
        tmp = zmax;
        zmax = zmin*2;
        estA = Methods.julius(data,net,zmax,logical(1));
    end
    % refine upper bound
    while zmax-zmin > tol
        i = (zmax + zmin) * 0.5;
        estA = Methods.julius(data,net,i,logical(1));
        if nnz(estA) == 0
            zmax = i;
        else
            zmin = i;
        end
    end

    % Convert to interval.
    delta = zetaRange(2)-zetaRange(1);
    zetavec = zetavec*delta + zetaRange(1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parse data to use
P = data.P;
Y = responce(data,net);

[cvY,cvP] = cov(data);

if ~isempty(data.cvY)
    cvY = data.cvY;
end

if ~isempty(data.cvP)
    cvP = data.cvP;
end

%% run

tol = 10^-6;
[nGenes, nExp] = size(data.P);

if ~exist('initA','var')
    initA = -eye(nGenes);
end

for i=1:length(zetavec)
    if size(initA,3) == length(zetavec)
        Ai = initA(:,:,i);
    else % we use one
        Ai = initA(:,:,1);
    end
    zeta = zetavec(i);
    [H, D, H] = svd(Ai * cvY * Ai' + cvP); % + e*eye(size(Ai)));
    dD = diag(D);
    % fprintf('sv, min = %g, max = %g\n',dD(end),dD(1))
    % Correction = H * sqrt(diag( 1./dD ));
    Correction = eye(size(Ai));

    cvx_begin quiet
    variables Ahat(nGenes,nGenes) eta(nGenes,nExp)
    minimize norm(eta'*Correction,'fro')+zeta*sum(sum(abs(Ahat)))
    subject to
    Ahat*Y+P == eta
    cvx_end
    
    S = (abs(Ahat) < tol);
    Ahat(S) = 0;
    Aest(:,:,i) = full(double(Ahat));
end



varargout{1} = Aest;
