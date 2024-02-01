function varargout = glmpoi(varargin)
% This is a wrapper function for the fitclinear:
%   logistic learner with lasso regularization
% function estA = lassolog(data,zetavec[net,alpha,rawZeta])
%
%   Input Arguments: lassolog(data,zetavec[net,alpha,rawZeta])
%   ================
%   data:    datastruct.Dataset
%   net:     datastruct.Network
%   zetavec: method parameter for tuning the network fitting.
%   alpha:   The elasticnet mixing parameter, with 0 < alpha <= 1. (default = 1, Lasso)
%            Currently alpha < 0.01 is not reliable, unless you
%            supply your own zeta sequence. zeta needs to be set first
%   rawZeta: logical to determine if the zeta values should be
%            converted.  default = false
%
%   Output Arguments: estA
%   =================
%   estA: the estimated networks as a 3d array.
%   fit: The glmnet algorithm output data, outputs the last fit if
%        a range of sparsity penalties are supplied.
%   zetaRange: the sparsity range used when a scaled zeta are used (default).
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse input arguments  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawZeta = 0;
zetavec = [];
net = [];
%alpha = 1;
for i = 1:nargin
    if isa(varargin{i},'datastruct.Dataset')
        data = varargin{i};
    elseif isa(varargin{i},'datastruct.Network')
        net = varargin{i};
    elseif isa(varargin{i},'logical')
        rawZeta = varargin{i};
    else
        if isempty(zetavec)
            zetavec = varargin{i};
        %else
            %alpha = varargin{i};
        end
    end
end

if ~exist('data','var')
    error('needs a data set')
end

%% Determine how to handle zeta %%

if ~rawZeta
    zetaRange = [];
    tol = 1e-6;
    zmax = 1;
    zmin = 0;
    % find zero network
    estA = Methods.glmpoi(data,net,zmax,true);
    while nnz(estA) > 0
        %tmp = zmax;
        zmax = zmax*2;
        estA = Methods.glmpoi(data,net,zmax,true);
    end
    % refine
    while zmax-zmin > tol
        i = (zmax + zmin) * 0.5;
        estA = Methods.glmpoi(data,net,i,true);
        if nnz(estA) == 0
            zmax = i;
        else
            zmin = i;
        end
    end

    zetaRange(1) = 0;
    zetaRange(2) = zmax;
    varargout{3} = zetaRange;
    % Convert to interval.
    delta = zetaRange(2) - zetaRange(1);
    zetavec = zetavec*delta + zetaRange(1);
end

%% Run
Afit = zeros(size(data.P,1),size(data.P,1),length(zetavec));
for i = 1:size(data.P,1)
   [Beta,fitinfo] = lassoglm(data.Y', -data.P(i,:)','poisson','Lambda',zetavec);
   Afit(i,:,:) = Beta(:,:);
end


varargout{1} = Afit;


if ~rawZeta
    varargout{2} = zetaRange;
    varargout{3} = fitinfo;
else
    varargout{2} = fitinfo;
end