function varargout = neunetcv(varargin)
%
% This is a wrapper function for the neural network algorihtm:
%
% function estA = nnet(data,zetavec,rawZeta)
%
%   Input Arguments: nnet(data,zetavec,rawZeta)
%   ================
%   data:    datastruct.Dataset
%   zetavec: method parameter for tuning the network fitting.
%   rawZeta: logical to determine if the zeta values should be
%            converted.  default = false
%
%   Output Arguments: estA
%   =================
%   estA: the estimated networks as a 3d array.
%   zetaRange: the sparsity range used when a scaled zeta are used (default).
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse input arguments  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawZeta = 0;
zetavec = [];

for i = 1:nargin
    if isa(varargin{i},'datastruct.Dataset')
        data = varargin{i};
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



%% Run
Afit = zeros(size(data.P,1),size(data.P,1));
    for i = 1:size(data.P,1)
        lambda=zetavec;
        cvloss=zeros(length(lambda),1);
        for j = 1:length(lambda)
            cvMdl = fitcnet(data.Y', data.P(i,:),"Activations","sigmoid","LayerSizes",1,Lambda=lambda(j), CVPartition=cvpartition(size(data.Y,2),"KFold",10));
            cvloss(j) = kfoldLoss(cvMdl);
        end
        [~, indmn] = min(cvloss);
        mdl = fitcnet(data.Y', data.P(i,:),"Activations","sigmoid","LayerSizes",size(data.Y,1),Lambda=zetavec(indmn));
        weis = mdl.LayerWeights{end};
        Afit(i,:) = weis(1,:);
    end

%% handle zeta %%

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
else
end
