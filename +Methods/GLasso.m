function varargout = GLasso(varargin)
% function estA = GLasso(data,net,zetavec,rawZeta)
%
%   Input Arguments: GLasso(data,zetavec[,rawZeta,net])
%   ================
%   data:    datastruct.Dataset
%   zetavec: method parameter for tuning the network fitting. (optional)
%   rawZeta: logical to determine if the zeta values should be
%            converted.  default = false
%   net:     datastruct.Network
%
%   Output Arguments: estA
%   =================
%   estA: the estimated networks as a 3d array.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse input arguments %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawZeta = false;

net = [];
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
    error('%s needs a data set',mfilename)
end

%% Determine how to handle zeta %%

if ~rawZeta
    zetaRange = [];
    tol = 1e-6;
    zmax = 1;
    zmin = 0;
    % find zero network
    estA = Methods.GLasso(data,net,zmax,logical(1));
    while nnz(estA) > 0 
        tmp = zmax;
        zmax = zmin*2;
        estA = Methods.GLasso(data,net,zmax,logical(1));
    end
    % refine
    while zmax-zmin > tol
        i = (zmax + zmin) * 0.5;
        estA = Methods.GLasso(data,net,i,logical(1));
        if nnz(estA) == 0
            zmax = i;
        else
            zmin = i;
        end
    end
    
    zetaRange(1) = 0;
    zetaRange(2) = zmax;
    varargout{2} = zetaRange;
    % Convert to interval.
    delta = zetaRange(2)-zetaRange(1);
    zetavec = zetavec*delta + zetaRange(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parse data to use
nY = response(data,net);

%% Run method with input paramters
Y = nY';

for i=1:length(zetavec)
    estA(:,:,i) = GraphicalLasso(Y, zetavec(i));
end

varargout{1} = estA;
