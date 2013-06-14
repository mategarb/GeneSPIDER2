function varargout = Glmnet(varargin)
% function estA = Glmnet(data,require,zetavec,rawZeta)
% evalGLasso will evaluate the method glasso and produce an
% estimated network matrix.
%
%   Input Arguments: Glmnet(data,net,zetavec,rawZeta)
%   ================
%   data:    GeneSpider.Dataset
%   net:     GeneSpider.Network
%   zetavec: method parameter for tuning the network fitting. (optinal)
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
    estA = Methods.Glmnet(data,net,zmax,logical(1));
    while nnz(estA) > 0 
        tmp = zmax;
        zmax = zmin*2;
        estA = Methods.Glmnet(data,net,zmax,logical(1));
    end
    % refine
    while zmax-zmin > tol
        i = (zmax + zmin) * 0.5;
        estA = Methods.Glmnet(data,net,i,logical(1));
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

%% Run
for i = 1:size(data.P,1)
    % fit = glmnet(nY',-P(i,:)','gaussian',glmnetSet(struct('lambda',zetavec,'standardize',false)));
    fit = glmnet(responce(data,net)',-data.P(i,:)','gaussian',glmnetSet(struct('lambda',zetavec)));
    Afit(i,:,:) = fit.beta(:,:);
end

Afit(:, :, :) = Afit(:, :, end:-1:1); % Glmnet reverses the order. Need to undo.

for i=1:size(Afit,3)
    Afit(:,:,i);
end

varargout{1} = Afit;

varargout{3} = fit;