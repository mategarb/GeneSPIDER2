function varargout = RNI(varargin)
% function estA = RNI(data <,net,alpha,regpath>)
%
%   Input Arguments: RNI(data <,net,alpha,zetavec,regpath>)
%   ================
%   data:    datastruct.Dataset
%   net:     datastruct.Network <optional>
%   alpha:   confidense level of inference, can only be a single digit \in [0,1[. default 0.01
%   zetavec: should the confidence values be truncated at some specified values.
%            Will be ignored if regpath = 'full'
%   regpath  {'input','full'}  string to determine if we should try to create a zetavec
%            for the complete regularization path dependant on method, default 'input'.
%            Where 'input' is a zetavec determined by the user
%            and 'full' gives the best estimate of the complete regularization path
%
%   Output Arguments: confA, zetavec
%   =================
%   confA: the confident networks as a 3d array.
%   zetavec: the complete regularization path
%   zetaRange: will return the raw zeta range scale factors
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse input arguments %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawZeta = 0;
regpath = 'input';
alpha = 0.01;
net = [];
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
        if length(varargin{i}) > 1
            zetavec = varargin{i};
        else
            alpha = varargin{i};
        end
    end
end

if ~exist('data')
    error('needs a data set')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = data.lambda;
if numel(lambda) == 2
    o = ones(size(data.P));
    [gamma, ps, Nrps] = RInorm(response(data,net)',data.P',diag(lambda(1:length(lambda)/2))*o',diag(lambda(length(lambda)/2+1:end))*o'+eps,alpha);
else
    o = ones(size(data.P),2);
    [gamma, ps, Nrps] = RInorm(response(data,net)',data.P',(diag(lambda(1:length(lambda)/2))'*o)',(diag(lambda(length(lambda)/2+1:end))'*o)'+eps,alpha);
end

if strcmpi(regpath,'full')
    zetavec = abs(gamma(:)');
    zetavec = unique(zetavec);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine how to handle zeta %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~rawZeta & strcmpi(regpath,'input')
    zetaRange = [];
    estA = gamma;
    zetaRange(1) = min(abs(estA(estA~=0)))-eps;
    zetaRange(2) = max(abs(estA(estA~=0)))+10*eps;

    % Convert to interval.
    delta = zetaRange(2)-zetaRange(1);
    zetavec = zetavec*delta + zetaRange(1);

elseif ~rawZeta & strcmpi(regpath,'full')
    zetaRange(2) = max(zetavec);
    zetaRange(1) = min(zetavec);
    delta = zetaRange(2)-zetaRange(1);
elseif rawZeta & strcmpi(regpath,'full')
    zetaRange(2) = 1;
    zetaRange(1) = 0;
    delta = 1;
end

for i=1:length(zetavec)
    temp = find(abs(gamma) <= zetavec(i));
    Atmp = gamma;
    Atmp(temp) = 0;

    estA(:,:,i) = Atmp;
end

varargout{1} = estA;

if nargout > 1 & strcmpi(regpath,'input')
    varargout{2} = zetavec;
    varargout{3} = zetaRange;
elseif nargout > 1 & strcmpi(regpath,'full')
    zetavec = (zetavec-zetaRange(1))/delta;
    varargout{2} = zetavec;
    varargout{3} = zetaRange;
end

return
