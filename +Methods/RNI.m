function varargout = RNI(varargin)
% function estA = RNI(data <,net>)
%
%   Input Arguments: RNI(data <,net>)
%   ================
%   data:    GeneSpider.Dataset
%   net:     GeneSpider.Network <optional>
%   alpha:   confidense level of inference.
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
        alpha = varargin{i}; % Obsolete
    end
end

if ~exist('data')
    error('needs a data set')
end

if ~exist('alpha')
    alpha = 0.01;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(alpha)
    lambda = data.lambda;
    if numel(lambda) == 2
        o = ones(size(data.P));
        [gamma, ps, Nrps] = tools.RInorm(response(data,net)',data.P',diag(lambda(1:length(lambda)/2))*o',diag(lambda(length(lambda)/2+1:end))*o'+eps,alpha(i));
    else
        o = ones(size(data.P),2);
        [gamma, ps, Nrps] = tools.RInorm(response(data,net)',data.P',(diag(lambda(1:length(lambda)/2))'*o)',(diag(lambda(length(lambda)/2+1:end))'*o)'+eps,alpha(i));
    end
    Gamma(:,:,i) = gamma;
    PS(:,:,i) = double(ps);
    NrPS(i) = Nrps;
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
