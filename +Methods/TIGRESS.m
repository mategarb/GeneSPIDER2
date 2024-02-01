function varargout = TIGRESS(varargin)
% function estA = tigress(data,R,L,alpha)
%
%   Input Arguments: tigress(data,R,L,alpha)
%   ================
%   data:    datastruct.Dataset
%   R:	     number of resamplings that should be used to run stability selection
%   	     (default=1000). Note that R can be a vector of increasing values. In
%   	     this case, TIGRESS will return frequencies for each of these values.
%   L:       number of LARS steps that should be considered
%   alpha:   randomization level. alpha is a scalar st 0<alpha<=1.
%	     If alpha=1, no randomization is used (default=.2)
%%
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
        dat = varargin{i};
    elseif isa(varargin{i},'datastruct.Network')
        net = varargin{i};
    elseif isa(varargin{i},'char')
	% zetavec = varargin{i};
	regpath = varargin{i};
    elseif isa(varargin{i},'logical')
        rawZeta = varargin{i};
    % elseif isa(varargin{i},'integer')
    %     restk = varargin{i};
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

data = struct('expdata', (dat.Y)');
% data = struct('expdata', (dat.Y));
tfindices = double(1:size(dat.Y,1));
% tfindices = double(randi([1 100],1,100));
% freq = tigress(data,'R', 4000, 'L',5, 'alpha', 0.2); % R 10000; L 2; a 0.4

% [edges scores freq] = tigress_full(data.expdata)


freq = tigress(data,'R', 1000, 'alpha', 0.4, 'L',2);
scores = score_edges(freq);
%scores = score_edges(freq, 'method', 'area', 'L', 3);
edges = predict_network(scores, tfindices);
%edges = predict_network(scores, tfindices, 'cutoff', 15);
% edges = predict_network(scores);
Network = zeros(size(dat.Y,1), size(dat.Y,1));
for i = 1:size(edges,1)
    Network(edges(i,2), edges(i,1)) = edges(i,3);
end


if ~exist('zetavec','var') & strcmpi(regpath,'input')
    estA = Network;
    %varargout{1} = estA;
    return
else
    Als = Network;
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
    elseif strcmpi(regpath, 'CO')
        zetavec = sort(abs(Als(:)'));
        zetavec = unique(zetavec);
        z = randperm(length(zetavec), 18);
        zetavec = [min(zetavec) sort(zetavec(z)) max(zetavec)];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(zetavec)
    temp = abs(Als) <= zetavec(i);
    Atmp = Als;
    Atmp(temp) = 0;
    % estA{i} = Atmp;
    estA(:,:,i) = Atmp;
    % estA = Atmp;
    % M{i} = analyse.CompareModels((estA), (Net.A));

end

% varargout{1} = M;

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
