function varargout = lasso(varargin)
% function estA = lasso(data,net,zetavec[,alpha,rawZeta])
%
%   Input Arguments: lasso(data,net,zetavec[,alpha,rawZeta])
%   ================
%   data:    datastruct.Dataset
%   net:     datastruct.Network
%   zetavec: method parameter for tuning the network fitting.
%   rawZeta: logical to determine if the zeta values should be
%            converted.  default = false
%
%   Output Arguments: estA
%   =================
%   estA: the estimated networks as a 3d array.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse input arguments %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawZeta = 0;
zetavec = [];
for i=1:nargin
    if isa(varargin{i},'datastruct.Dataset')
        data = varargin{i};
    elseif isa(varargin{i},'datastruct.Network')
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
if ~exist('net')
    X = (data.Y+data.E)';
else
    X = response(data,net)';
end

%% Determine how to handle zeta %%
% X = normalize(response(data,net)')
if ~rawZeta
    zmax = 0;
    for i = 1:size(data.P,1)
        % y = center(-data.P(i,:)');
        y = -data.P(i,:)';
        [beta, info] = lasso(X,y,0,true); % calculates leaste squares solution
        if sum(abs(beta(:,end))) > zmax
            zmax = sum(abs(beta(:,end)));
        end
    end

    zetaRange = [0, zmax];

    % Convert to interval.
    delta = zetaRange(2)-zetaRange(1);
    zetavec = zetavec*delta + zetaRange(1);

    varargout{3} = zetaRange;
else
    varargout{3} = [min(zetavec), max(zetavec)];
end

%% Run
for j=1:length(zetavec)
    for i = 1:size(data.P,1)
        % y = center(-data.P(i,:)')
        y = -data.P(i,:)';
        [beta, info] = lasso(X,y,zetavec(j),false);
        Afit(i,:,j) = beta';
    end
end
Afit(:, :, :) = Afit(:, :, end:-1:1); % lasso reverses the order. Need to undo.

varargout{1} = Afit;

varargout{2} = info;