function varargout = RNI(varargin)
% function estA = RNI(data <,net>)
%
%   Input Arguments: RNI(data <,net>)
%   ================
%   data:    GeneSpider.Dataset
%   net:     GeneSpider.Network <optional>
%
%   Output Arguments: confA
%   =================
%   confA: the confident networks as a 3d array.
%   

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
        rawZeta = varargin{i}; % Obsolete
    else
        zetavec = varargin{i}; % Obsolete
    end
end

if ~exist('data')
    error('needs a data set')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine how to handle zeta %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if ~rawZeta
%     zetaRange = [];
%     estA = -data.P*pinv(response(data,net));
%     zetaRange(1) = min(abs(estA(estA~=0)))-eps;
%     zetaRange(2) = max(abs(estA(estA~=0)))+10*eps;
    
%     % Convert to interval.
%     delta = zetaRange(2)-zetaRange(1);
%     zetavec = zetavec*delta + zetaRange(1);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = data.lambda;
if numel(lambda) == 2
    o = ones(size(data.P));
    [Gamma, PS, NrPS] = tools.RInorm(response(data,net)',data.P',diag(lambda(1:length(lambda)/2))*o',diag(lambda(length(lambda)/2+1:end))*o'+eps,data.alpha);
else
    o = ones(size(data.P),2);
    [Gamma, PS, NrPS] = tools.RInorm(response(data,net)',data.P',(diag(lambda(1:length(lambda)/2))'*o)',(diag(lambda(length(lambda)/2+1:end))'*o)'+eps,data.alpha);
end

if nargout == 0,
    disp(Gamma);
elseif nargout == 1;
    varargout{1} = PS; 
elseif nargout == 2;
    varargout{1} = PS; 
    varargout{2} = Gamma; 
elseif nargout == 3;
    varargout{1} = PS; 
    varargout{2} = Gamma; 
    varargout{3} = NrPS; 
end
return
