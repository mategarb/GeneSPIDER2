function varargout = tlsco(varargin)
% function estA = tlsco(data,require,zetavec,rawZeta)
%
%   Input Arguments: tlsco(data,net,zetavec,rawZeta)
%   ================
%   data:    datastruct.Dataset
%   net:     datastruct.Network
%   zetavec: vector containing zeta values. In the context of tlsco
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
    hatTheta = Methods.tls(response(data,net)',data.P');
    estA = hatTheta';
    varargout{1} = estA;
else
    Atls = Methods.tls(response(data,net)',data.P');
end

if strcmpi(regpath,'full')
    zetavec = abs(Atls(:)');
    zetavec = unique(zetavec);
    % zetavec = zetavec(zetavec~=0)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine how to handle zeta %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~rawZeta & strcmpi(regpath,'input')
    zetaRange = [];
    hatTheta = Methods.tls(response(data,net)',-data.P');
    % Atls = hatTheta';
    estA = hatTheta'; %Atls
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

hatTheta = Methods.tls(response(data,net)',-data.P');
Atls = hatTheta';
for i=1:length(zetavec)
    temp = find(abs(Atls) <= zetavec(i));
    Atmp = Atls;
    Atmp(temp) = 0;

    estA(:,:,i) = Atmp;
    zetaRange(1) = min(abs(estA(estA~=0)))-eps;
    zetaRange(2) = max(abs(estA(estA~=0)))+10*eps;

end

varargout{1} = estA;


varargout{2} = zetaRange;

return
