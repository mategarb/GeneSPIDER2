function varargout = ridgeco(varargin)
% This function performes a ridge regression fit using the glmnet algorihtm:
% (Glmnet for Matlab (2013) Qian, J., Hastie, T., Friedman, J., Tibshirani, R. and Simon, N.
% http://www.stanford.edu/~hastie/glmnet_matlab/)
%
% function estA = ridgeco(data,zetavec[net,alpha,rawZeta])
%
%   Input Arguments: Glmnet(data,zetavec[net,alpha,rawZeta])
%   ================
%   data:    datastruct.Dataset
%   net:     datastruct.Network
%   zetavec: method parameter for tuning the network fitting.
%   rawZeta: logical to determine if the zeta values should be
%            converted.  default = false
%   regpath  {'input','full'}  string to determine if we should try to create a zetavec
%            for the complete regularization path dependant on method, default 'input'.
%            Where 'input' is a zetavec determined by the user
%            and 'full' lets the glmnet algorithm determine the relevant regularization penalties.
%            If the SNR is high this might lead to few peanalty steps as the first small peanalty
%            removes a majority of the interactions.
%
%   Output Arguments: estA
%   =================
%   estA: the estimated networks as a 3d array.
%   L2_pen: the penelty value with the highest Rsqr value giving
%   the full network used to apply the cutoff to.
%   zetavec: the complete regularization path
%   zetaRange: the sparsity range used when a scaled zeta are used (default).
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse input arguments  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawZeta = 0;
zetavec = [];
net = [];
alpha = 0;
regpath = 'input';
for i=1:nargin
    if isa(varargin{i},'datastruct.Dataset')
        data = varargin{i};
    elseif isa(varargin{i},'datastruct.Network')
        net = varargin{i};
    elseif isa(varargin{i},'logical')
        rawZeta = varargin{i};
    elseif isa(varargin{i},'char')
        regpath = varargin{i};
    else
        zetavec = varargin{i};
    end
end

if ~exist('data')
    error('needs a data set')
end

% perform the fit between pertubation and data 
% due to how glmnet works this have to be done 1 pertubation row at
% a time
% we run for glmnets default 100 lambdas and then we pick the
% highest Rsqr of this as the network to cut from. 
for i = 1:size(data.P,1)
    fit = glmnet(data.Y',-data.P(i,:)','gaussian',glmnetSet(struct('alpha',alpha)));
    Afit(i,:,:) = fit.beta(:,:);
end
% Glmnet reverses the order so reverse it again to get expected pattern
Als(:, :, :) = Afit(:, :, end:-1:1); 

TSS = sum((data.Y-mean(data.Y)).^2);

for i = 1:size(Als,3)
    yfit = -pinv(Als(:,:,i))*data.P;
    RSS = sum((data.Y-yfit).^2);
    Rsquared(i) = 1 - RSS/TSS;
end
% select the highest scoring dataset 
[~,select] = max(Rsquared);
L2_pen = flip(fit.lambda);
L2_pen = L2_pen(select);
% and set this to the new Als
Als = Als(:,:,select);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine how to handle zeta for cutoff %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(regpath,'full')
    zetavec = abs(Als(:)');
    zetavec = unique(zetavec);
end
if ~rawZeta & strcmpi(regpath,'input')
    zetaRange = [];
    estA = Als;
    zetaRange(1) = min(abs(estA(estA~=0)))-eps;
    zetaRange(2) = max(abs(estA(estA~=0)))+10*eps;

    % Convert to interval.
    delta = zetaRange(2)-zetaRange(1);
    zetavec = zetavec*delta + zetaRange(1);
elseif ~rawZeta & strcmpi(regpath,'full')
    zetavec=zetavec;
    zetaRange(2) = max(zetavec);
    zetaRange(1) = min(zetavec);
    delta = zetaRange(2)-zetaRange(1);
elseif rawZeta & strcmpi(regpath,'full')
    zetavec=zetavec;
    zetaRange(2) = 1;
    zetaRange(1) = 0;
    delta = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop over each sparsity and each matrix from the fit 
% and remove values smaller than zeta
for i=1:length(zetavec)
    temp = find(abs(Als) <= zetavec(i));
    Atmp = Als;
    Atmp(temp) = 0;

    estA(:,:,i) = Atmp;
end

% finally capture all output arguments 
varargout{1} = estA;

if nargout > 1 & strcmpi(regpath,'input')
    varargout{2} = L2_pen;
    varargout{3} = zetavec;
    varargout{4} = zetaRange;
elseif nargout > 1 & strcmpi(regpath,'full')
    zetavec = (zetavec-zetaRange(1))/delta;
    varargout{2} = L2_pen;
    varargout{3} = zetavec;
    varargout{4} = zetaRange;
end

return
