function varargout = lsco(varargin)
% function estA = lsco(data,require,zetavec,rawZeta)
%
%   Input Arguments: lsco(data,net,zetavec,rawZeta)
%   ================
%   data:    GeneSpider.Dataset
%   net:     GeneSpider.Network
%   zetavec: vector containing zeta values. In the context of lsco
%            this is a threshold value of the range from min to max element
%   rawZeta: logical to determine if the zeta values should be
%            converted.  default 0
%
%   Output Arguments: estA
%   =================
%   estA: the estimated networks as a 3d array.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse input arguments %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawZeta = 0;

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

%% Run
if ~exist('zetavec','var')
    estA = -data.P*pinv(response(data,net));
    varargout{1} = estA;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine how to handle zeta %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~rawZeta
    zetaRange = [];
    estA = -data.P*pinv(response(data,net));
    zetaRange(1) = min(abs(estA(estA~=0)))-eps;
    zetaRange(2) = max(abs(estA(estA~=0)))+10*eps;
    
    % Convert to interval.
    delta = zetaRange(2)-zetaRange(1);
    zetavec = zetavec*delta + zetaRange(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Als = -data.P*pinv(response(data,net));
for i=1:length(zetavec)
    temp = find(abs(Als) <= zetavec(i));
    Atmp = Als;
    Atmp(temp) = 0;

    estA(:,:,i) = Atmp;
end

varargout{1} = estA;
return
