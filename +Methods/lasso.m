function varargout = lasso(varargin)
% This is a wrapper function for the glmnet algorihtm:
% (Glmnet for Matlab (2013) Qian, J., Hastie, T., Friedman, J., Tibshirani, R. and Simon, N.
% http://www.stanford.edu/~hastie/glmnet_matlab/)
%
% function estA = Glmnet(data,zetavec[net,alpha,rawZeta])
%
%   Input Arguments: Glmnet(data,zetavec[net,alpha,rawZeta])
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
alpha = 1;
for i=1:nargin
    if isa(varargin{i},'datastruct.Dataset')
        data = varargin{i};
    elseif isa(varargin{i},'datastruct.Network')
        net = varargin{i};
    elseif isa(varargin{i},'logical')
        rawZeta = varargin{i};
    else
        if isempty(zetavec)
            zetavec = varargin{i};
        else
            alpha = varargin{i};
        end
    end
end

if ~exist('data')
    error('needs a data set')
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
        zmax = zmax*2;
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
    varargout{3} = zetaRange;
    % Convert to interval.
    delta = zetaRange(2)-zetaRange(1);
    zetavec = zetavec*delta + zetaRange(1);
end

%% Run
for i = 1:size(data.P,1)
    % fit = glmnet(nY',-P(i,:)','gaussian',glmnetSet(struct('lambda',zetavec,'standardize',false)));
    fit = glmnet(response(data,net)',-data.P(i,:)','gaussian',glmnetSet(struct('lambda',zetavec,'alpha',alpha)));
    Afit(i,:,:) = fit.beta(:,:);
end
Afit(:, :, :) = Afit(:, :, end:-1:1); % Glmnet reverses the order. Need to undo.

% for i=1:size(Afit,3)
%     Afit(:,:,i);
% end

varargout{1} = Afit;


if ~rawZeta
    varargout{2} = zetaRange;
    varargout{3} = fit;
else
    varargout{2} = fit;
end
