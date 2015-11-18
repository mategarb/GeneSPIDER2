function varargout = Bolsco(varargin)
% Bootstrap sampling applied for least squares with cut off
%
% function [estA,Asupport] = Bolsco(data,bootstraps,zetavec[,net,rawZeta,straps])
%
%   Input Arguments:
%   ================
%   data:    datastruct.Dataset
%   net:     datastruct.Network
%   zetavec: vector containing zeta values. In the context of lsco
%            this is a threshold value of the range from min to max element
%   bootstraps: number of bootstraps, (default = 100, will not accept = 1)
%   rawZeta: logical to determine if the zeta values should be
%            converted.  default 0
%   straps:  integer, number of bootstrap runs. must be > 1, default = 100
%
%   Output Arguments: estA
%   =================
%   estA: the estimated networks as a 3d array (100% bootstrap supported network, binary).
%   Asupport: the support after <straps> bootstraps.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse input arguments  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawZeta = 0;
tmpzetas = [];
net = [];
alpha = 1;
tmpstraps = 1;
for i=1:nargin
    if isa(varargin{i},'datastruct.Dataset')
        data = varargin{i};
    elseif isa(varargin{i},'datastruct.Network')
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

zR = []; % zeta range for all bootstraps
Apos = zeros(data.N,data.N,straps);
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
        zR(:,j) = zetaRange;
        varargout{3} = zR;

        % Convert to interval.
        delta = zetaRange(2)-zetaRange(1);
        zetavec = zetavec*delta + zetaRange(1);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Als = -bdata.P*pinv(response(bdata,net));
    for i=1:length(zetavec)
        temp = find(abs(Als) <= zetavec(i));
        Atmp = Als;
        Atmp(temp) = 0;
        %disp(nnz(Atmp))
        estA(:,:,i) = Atmp;
    end

    if isempty(Alogical)
        Alogical = double(logical(estA));
        Asign=(double(sign(estA))<0);
        tmp = estA > 0;
        Apos = double(tmp);
        tmp = estA < 0;
        Aneg = double(tmp);
        Aneg = abs(Aneg);
    else
        Alogical = Alogical + double(logical(estA));
        tmp = estA > 0;
        Apos = Apos + double(tmp);
        tmp = estA < 0;
        Aneg = Aneg + double(tmp);
    end

end

Alogical = Alogical/straps;
Afrac = Alogical;
Alogical(Alogical < 1) = 0;
Aposs_frac = Apos./(straps*Afrac);
Asign_frac = 2*Aposs_frac-1;

for i=1:size(Apos,3)
    tmp=cat(3,Apos,Aneg);
    Asign_Max(:,:,i)=max(tmp,[],3);
end
Asign_Max=(Asign_Max.*Afrac)/straps;

if nargout > 0
    varargout{1} = Alogical;
end

if nargout > 1
    varargout{2} = Afrac;
end

if nargout > 2 % Agnostic sign support
    varargout{3} = Asign_Max;
end

if nargout > 3 % Explicit sign support -1 for 100% negative
    varargout{4} = Asign_frac;
end

if nargout > 4 % # positive links
    varargout{5} = cat(3,Apos,Aneg);
end

return
