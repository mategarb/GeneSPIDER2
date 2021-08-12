function varargout = LSCON(varargin)
% function estA = lsco(data,net,zetavec,rawZeta,regpath)
%
%   Input Arguments: lsco(data,net,zetavec,rawZeta,regpath)
%   ================
%   data:    datastruct.Dataset
%   net:     datastruct.Network
%   zetavec: vector containing zeta values. In the context of lsco
%            this is a threshold value of the range from min to max element.
%            Will be ignored if regpath = 'full' is specified
%            If regpath='AZ' this should be an int with how many
%            steps should be in the zetavec that are generated
%   rawZeta: logical to determine if the zeta values should be
%            converted.  default 0
%   regpath  {'input','full','AZ'}  string to determine if we should try to create a zetavec
%            for the complete regularization path dependant on method, default 'input'.
%            Where 'input' is a zetavec determined by the user,
%            'full' gives the best estimate of the complete
%            regularization path
%            'AZ' draws a zetavec that spans the full network and
%            the emtpy network in equal steps with each step being
%            n links smaller than the last one based on the steps
%            determined by the zetavec
%
%   Output Arguments: estA, zetavec, zetaRange
%   =================
%   estA: the estimated networks as a 3d array.
%   zetavec: if regpath is set to 'full' will return the used normalized zetavec
%   zetaRange: will return the raw zeta range scale factors
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse input arguments %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawZeta = 0;
regpath = 'input';
net = [];
for i=1:nargin
    if isa(varargin{i},'datastruct.Dataset')
        data = varargin{i};
    elseif isa(varargin{i},'datastruct.Network')
        net = varargin{i};
    elseif isa(varargin{i},'char')
        regpath = varargin{i};
    elseif isa(varargin{i},'logical')
        rawZeta = varargin{i};
    else
        zetavec = varargin{i};
    end
end

if ~exist('data')
    error('needs a data set')
end

%% Run
if ~exist('zetavec','var') & strcmpi(regpath,'input')
    estA = -data.P*pinv(response(data,net));
    % run normalization 
    for i = 1:size(estA,2)
        estA(:,i) = estA(:,i)./sum(abs(estA(:,i)))*data.N;
    end
    varargout{1} = estA;
    return
else
    Als = -data.P*pinv(response(data,net));
    % run normalization 
    for i = 1:size(Als,2)
        Als(:,i) = Als(:,i)./sum(abs(Als(:,i)))*data.N;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine how to handle zeta %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(regpath,"AZ")
    % test if zetavec is a array, if yes get the length
    if length(zetavec) > 1
        zetavec = length(zetavec);
    end
    % then draw a set of values that give a stepwise sparsity
    % between full and empty
    zetavec = gsUtilities.adaptive_sparsity(Als,zetavec);
else
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
        zetaRange(2) = max(zetavec);
        zetaRange(1) = min(zetavec);
        delta = zetaRange(2)-zetaRange(1);
    elseif rawZeta & strcmpi(regpath,'full')
        zetaRange(2) = 1;
        zetaRange(1) = 0;
        delta = 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(zetavec)
    temp = find(abs(Als) <= zetavec(i));
    Atmp = Als;
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
