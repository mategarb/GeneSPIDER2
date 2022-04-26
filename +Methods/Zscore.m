function varargout = Zscore(varargin)
%
%   Input Arguments: Zscore(data, zetavec)
%   ================
%   data:    datastruct.Dataset
%   zetavec: 'CO' if one wants to get 20 network with different sparstiy levels,
%            'full' if one wants to get all possible sparsities
%   Output Arguments: estA, zetavec, zetaRange
%   =================
%   estA: the estimated networks as a 3d array.
%   zetavec: if regpath is set to 'full' or 'CO' will return the used normalized zetavec
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
        dat = varargin{i};
    elseif isa(varargin{i},'datastruct.Network')
        net = varargin{i};
    elseif isa(varargin{i},'char')
	     regpath = varargin{i};
    elseif isa(varargin{i},'logical')
        rawZeta = varargin{i};
     else
        if length(varargin{i}) > 1 | varargin{i} > 1
            zetavec = varargin{i};
        else
            alpha = varargin{i};
        end
    end
end

if ~exist('dat')
    error('needs a data set')
end

%% Run

Zsc = zeros(dat.N,dat.N);
for z = 1:size(dat.Y,1)
  loc = find(dat.P(z,:));
  d = dat.Y(:,loc); ave = mean(d,2);
  % ind = 1:dat.N;
  for w = 1:length(ave)
    % ig = find(ind ~= ind(w));
    Zsc(w,z) = (ave(w) - mean(dat.Y(w,:)))/std(dat.Y(w,:));
  end
end

if ~exist('zetavec','var') & strcmpi(regpath,'input')
    estA = Zsc;
    return
else
    Als = Zsc;
end

if strcmpi(regpath,'AZ')
    % test if zetavec is a array, if yes get the length
    if length(zetavec) > 1
        zetavec = length(zetavec);
    end
    % then draw a set of values that give a stepwise sparsity
    % between full and empty
    zetavec = gsUtilities.adaptive_sparsity(Als,zetavec);
elseif strcmpi(regpath,'full')
    zetavec = sort(abs(Als(:)'));
    zetavec = unique(zetavec);
elseif strcmpi(regpath, 'CO')
    zetavec = sort(abs(Als(:)'));
    zetavec = unique(zetavec);
    z = randperm(length(zetavec), 18);
    zetavec = [min(zetavec) sort(zetavec(z)) max(zetavec)];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine how to handle zeta %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
