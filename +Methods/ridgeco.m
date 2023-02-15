function varargout = ridgeco(varargin)
%
% This is a wrapper function for ridge from matlab toolbox
% function estA = mridgeco(data,net,zetavec,rawZeta,regpath)
%
%   Input Arguments: mridgeco(data,net,zetavec,rawZeta,regpath)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse input arguments  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawZeta = 0;
zetavec = [];
net = [];
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

if ~exist('data','var')
    error('Data set is missing')
end

%% run ridge

Afit = zeros(size(data.P,1),size(data.P,1),length(zetavec));
    for i = 1:size(data.P,1)
        
        b = ridge(data.P(i,:)',response(data,net)', zetavec, 1);
        Afit(i,:,:) = b(:,:);

        % some other ways to run ridge
        % [b, stats] = lasso(response(data,net)', -data.P(i,:)', 'Lambda', zetavec,'Alpha',eps(1));
        % [b, stats] = fitrlinear(response(data,net)', -data.P(i,:)','Learner','leastsquares','Lambda', zetavec,'Regularization','ridge');

    end

TSS = sum((data.Y-mean(data.Y)).^2);
Rsquared = zeros(size(Afit,3),1);
for i = 1:size(Afit,3)
    yfit = -(pinv(Afit(:,:,i))*data.P);
    RSS = sum((data.Y-yfit).^2);
    Rsquared(i) = 1 - RSS/TSS;
end
% select the highest scoring dataset 
[~,select] = max(Rsquared);
L2_pen = flip(zetavec);
L2_pen = L2_pen(select);
% and set this to the new Als
Als = Afit(:,:,select);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine how to handle zeta for cutoff %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(regpath,'full')
    zetavec = abs(Als(:)');
    zetavec = unique(zetavec);
end
if ~rawZeta && strcmpi(regpath,'input')
    zetaRange = [];
    estA = Als;
    zetaRange(1) = min(abs(estA(estA~=0)))-eps;
    zetaRange(2) = max(abs(estA(estA~=0)))+10*eps;

    % Convert to interval.
    delta = zetaRange(2)-zetaRange(1);
    zetavec = zetavec*delta + zetaRange(1);
elseif ~rawZeta && strcmpi(regpath,'full')
    zetaRange(2) = max(zetavec);
    zetaRange(1) = min(zetavec);
    delta = zetaRange(2)-zetaRange(1);
elseif rawZeta && strcmpi(regpath,'full')
    zetaRange(2) = 1;
    zetaRange(1) = 0;
    delta = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop over each sparsity and each matrix from the fit 
% and remove values smaller than zeta
for i = 1:length(zetavec)
    Atmp = Als;
    Atmp(abs(Als) <= zetavec(i)) = 0;
    estA(:,:,i) = Atmp;
end

% finally capture all output arguments 
varargout{1} = estA;

if nargout > 1 && strcmpi(regpath,'input')

    varargout{2} = L2_pen;
    varargout{3} = zetavec;
    varargout{4} = zetaRange;

elseif nargout > 1 && strcmpi(regpath,'full')

    zetavec = (zetavec-zetaRange(1))/delta;
    varargout{2} = L2_pen;
    varargout{3} = zetavec;
    varargout{4} = zetaRange;

end

return
