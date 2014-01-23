function varargout = ccd(varargin)
% function estA = ccd(data,net,zetavec,rawZeta)
%
%   Input Arguments: ccd(data,net,zetavec,rawZeta)
%   ================
%   data:    GeneSpider.Dataset
%   net:     GeneSpider.Network
%   zetavec: method parameter for tuning the network fitting. (optional)
%   rawZeta: logical to determine if the zeta values should be
%            converted.  default = false
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

%% Determine how to handle zeta %%

if ~rawZeta
    zetaRange = [];
    tol = 1e-6;
    zmax = 1;
    zmin = 0;
    % find zero network
    estA = Methods.ccd(data,net,zmax,logical(1));
    while nnz(estA) > 0 
        tmp = zmax;
        zmax = zmin*2;
        estA = Methods.ccd(data,net,zmax,logical(1));
    end
    % refine
    while zmax-zmin > tol
        i = (zmax + zmin) * 0.5;
        estA = Methods.ccd(data,net,i,logical(1));
        if nnz(estA) == 0
            zmax = i;
        else
            zmin = i;
        end
    end
    
    zetaRange(1) = 0;
    zetaRange(2) = zmax;
    varargout{2} = zetaRange;
    % Convert to interval.
    delta = zetaRange(2)-zetaRange(1);
    zetavec = zetavec*delta + zetaRange(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parse data to use.
Xt = responce(data,net);
[m,n] = size(Xt'); % n variables, m equations
X = Xt' - repmat(mean(Xt',1),size(Xt',1),1);

if isfield(data,'P'),
    Ut = data.P;
    U = Ut' - repmat(mean(Ut',1),size(Ut',1),1);
    b = zeros(n,1); 
    for k = 1:n; b(k)=regress(X(:,k), U(:,k)); end;
else
    U = zeros(size(X));
    b = zeros(n,1);
end

%% Adjusting for U and finding zeta_max
Y = zeros(size(X));
for k = 1:n
  Y(:,k) = X(:,k) - b(k)*U(:,k);
end

qk = 1;

XtX = X' * X;

%% Run =============================================

for i=1:length(zetavec)
    C = sparse(zeros(n,n));
    ks = 1;
    zeta = zetavec(i);
    while ks <= n
        beta=ccd2(X,Y(:,ks),zeta,'XtX',XtX,'trace',1);   
        C(ks,:) = beta;
        ks = ks + 1;
    end
    
    Aest(:,:,i) = full(C);
end

varargout{1} = Aest;
