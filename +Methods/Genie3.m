function varargout = Genie3(varargin)

rawZeta = 0;
regpath='input';
net = [];

for i=1:nargin
    if isa(varargin{i},'datastruct.Dataset')
        expr_matrix = varargin{i};
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

if ~exist('expr_matrix')
    error('%s needs a data set', mfilename)
end

if ~exist('zetavec','var') & strcmpi(regpath, 'input')
    estA = genie3(expr_matrix.Y');
    ns = size(expr_matrix.Y);
    estA = estA'+eye(ns(1));
else
    Als = genie3(expr_matrix.Y');
    ns = size(expr_matrix.Y);
    Als = Als'+eye(ns(1));
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

for i=1:length(zetavec)
    temp = find(abs(Als) <= zetavec(i));
    Atmp = Als;
    Atmp(temp) = 0;
    estA(:,:,i) = Atmp';
end
varargout{1} = estA;

if nargout > 1 & strcmpi(regpath, 'input')
    varargout{2} = zetavec;
 %   varargout{3} = zetaRange;
elseif nargout > 1 & strcmpi(regpath, 'full')
    zetavec = (zetavec-zetaRange(1))/delta;
    varargout{2} = zetavec;
  %  varargout{3} = zetaRange;
end
return
