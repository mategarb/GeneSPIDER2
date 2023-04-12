function varargout = svm(varargin)

% This is a wrapper function for the support vector machine classifier:
%
% function estA = svm(data,zetavec[net,rawZeta])
%
%   Input Arguments: Glmnet(data,zetavec[net,rawZeta])
%   ================
%   data:    datastruct.Dataset
%   net:     datastruct.Network
%   zetavec: method parameter for tuning the network fitting.
%   rawZeta: logical to determine if the zeta values should be
%            converted.  default = false
%
%   Output Arguments: estA
%   =================
%   estA: the estimated networks as a 3d array.
%   fit: The glmnet algorithm output data, outputs the last fit if
%        a range of sparsity penalties are supplied.
%   zetaRange: the sparsity range used when a scaled zeta are used (default).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse input arguments  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rawZeta = 0;
zetavec = [];
%net = [];

for i = 1:nargin
    if isa(varargin{i},'datastruct.Dataset')
        data = varargin{i};
    elseif isa(varargin{i},'datastruct.Network')
        %net = varargin{i};
    elseif isa(varargin{i},'logical')
        rawZeta = varargin{i};
    else
        if isempty(zetavec)
            zetavec = varargin{i};
        end
    end
end

if ~exist('data','var')
    error('needs a data set')
end

%% Run svm

Afit = zeros(size(data.P,1),size(data.P,1));

    for i = 1:size(data.P,1)
        mdl = fitcsvm(data.Y', -data.P(i,:),'KernelFunction', 'linear'); % for low and moderate dimensions
        Afit(i,:) = mdl.Beta;
    end
%% Determine how to handle zeta %%

    zetaRange = [];
    estA = Afit;
    zetaRange(1) = min(abs(estA(estA~=0)))-eps;
    zetaRange(2) = max(abs(estA(estA~=0)))+10*eps;

    % Convert to interval.
    delta = zetaRange(2)-zetaRange(1);
    zetavec = zetavec*delta + zetaRange(1);

for i = 1:length(zetavec)
    Atmp = Afit;
    Atmp(abs(Afit) <= zetavec(i)) = 0;
    estA(:,:,i) = Atmp;
end

    varargout{1} = estA;

if ~rawZeta
    varargout{2} = zetaRange;


end
