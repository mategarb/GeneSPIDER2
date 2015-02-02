function varargout = Bolsco(varargin)
% Bootstrap sampling applied for least squares with cut off
%
% function [estA,Asupport] = Bolsco(data,bootstraps,zetavec{,net,rawZeta})
%
%   Input Arguments:
%   ================
%   data:    GeneSpider.Dataset
%   net:     GeneSpider.Network
%   zetavec: vector containing zeta values. In the context of lsco
%            this is a threshold value of the range from min to max element
%   bootstraps: number of bootstraps, (default = 100, will not accept = 1)
%   rawZeta: logical to determine if the zeta values should be
%            converted.  default 0
%
%   Output Arguments: estA
%   =================
%   estA: the estimated networks as a 3d array (100% bootstrap supported network, binary).
%   Asupport: the support after <straps> bootstraps.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse input arguments %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawZeta = 0;
tmpzetas = [];
net = [];
alpha = 1;
tmpstraps = 1;
for i=1:nargin
    if isa(varargin{i},'GeneSpider.Dataset')
        data = varargin{i};
    elseif isa(varargin{i},'GeneSpider.Network')
        net = varargin{i};
    elseif isa(varargin{i},'logical')
        rawZeta = varargin{i};
    elseif varargin{i}==floor(varargin{i}) & length(varargin{i}) == 1 & varargin{i} ~= 1 % is integer of length 1 and is not 1
        tmpstraps = varargin{i};
    else
        if isempty(tmpzetas)
            tmpzetas = varargin{i};
        else
            alpha = varargin{i};
        end
    end
end

if ~exist('data')
    error('needs a data set')
end

if tmpstraps == 1
    straps = 100;
else
    straps = tmpstraps;
end

%% Run
if ~exist('tmpzetas','var')
    estA = -data.P*pinv(response(data,net));
    varargout{1} = estA;
    return
end

Alogical = [];
for j = 1:straps
    zetavec = tmpzetas;
    bdata = bootstrap(data);

    reps = 1;
    while rank(bdata.P) < min(bdata.N,bdata.M)
        bdata = bootstrap(data);
        if ~mod(reps,10000)
            disp('number of booot tries')
            disp(reps)
        end
        reps = reps+1;
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

    if isempty(Alogical)
        Alogical = double(logical(estA));
    else
        Alogical = Alogical + double(logical(estA));
    end

end

Alogical = Alogical/straps;
Afrac = Alogical;
Alogical(Alogical < 1) = 0;
varargout{1} = Alogical;

if nargout > 1
    varargout{2} = Afrac;
end

return
