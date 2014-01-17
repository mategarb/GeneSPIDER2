function varargout = NIR(varargin)
% function estA = NIR(data,net,zetavec[,rawZeta])
% 
%   Input Arguments: NIR(data,net,zetavec[,rawZeta])
%   ================
%   data:    GeneSpider.Dataset
%   net:     GeneSpider.Network
%   zetavec: method parameter for tuning the network fitting.
%   rawZeta: logical to determine if the zeta values should be
%            converted.  default = false
%
%   Output Arguments: estA, sdA, covA, DIAG
%   =================
%   estA: the estimated networks as a 3d array.
%   sdA: standard deviation of estA (optional)
%   covA: covariance matrix of estA (optional)
%   DIAG: diagnostics on fit (optional)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse input arguments %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawZeta = 0;
zetavec = [];
net = [];
for i=1:nargin
    if isa(varargin{i},'GeneSpider.Dataset')
        data = varargin{i};
    elseif isa(varargin{i},'GeneSpider.Network')
        net = varargin{i};
    elseif isa(varargin{i},'logical')
        rawZeta = varargin{i};
    else
        zetavec = varargin{i};
    end
end

if ~exist('data')
    error('needs a data set')
end

%% Determine how to handle zeta %%
if ~rawZeta,

    zetaRange = [0 data.N];
    % Convert to interval.
    zetavec = 1-zetavec;
    delta = zetaRange(2)-zetaRange(1);
    zetavec = ceil(zetavec*delta + zetaRange(1));
    if length(zetavec) == 0
        error('NIR could not resolve the zeta values')
    end
    [zetavec,loc,orig] = unique(zetavec);
else
    [zetavec,loc,orig] = unique(zetavec);
end


%% Parse data to use
nY = response(data,net);
P = data.P;

[sdY,sdP] = std(data);

if ~isempty(data.sdY)
    sdY = data.sdY;
end
if ~isempty(data.sdP)
    sdP = data.sdP;
end

%% Run
estA = [];
for i=1:length(zetavec)
    [estA(:,:,i), sdA{i}, covA{i}, DIAG(i)] = nir(nY,sdY,P,sdP,zetavec(i));
end

varargout{1} = estA;
varargout{2} = sdA;
varargout{3} = covA;
varargout{4} = DIAG;
