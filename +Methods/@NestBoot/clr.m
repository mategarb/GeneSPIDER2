function varargout = clr(varargin)
% Bootstrap sampling applied for least squares with cut off
%
% function [estA,Asupport] = Bolsco(data,[bootstraps,zetavec,net,rawZeta,straps,regpath])
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
%   regpath  {'input','full'}  string to determine if we should try to create a zetavec
%            for the complete regularization path dependant on method, default 'input'.
%            Where 'input' is a zetavec determined by the user
%            and 'full' lets the glmnet algorithm determine the relevant regularization penalties.
%            If the SNR is high this might lead to few peanalty steps as the first small peanalty
%            removes a majority of the interactions.
%
%   Output Arguments: estA
%   =================
%   estA: the estimated networks as a 3d array (100% bootstrap supported network, binary).
%   Asupport: the support after <straps> bootstraps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse input arguments  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawZeta = 0;
zetavec = [];
net = [];
alpha = 1;
tmpstraps = 1;
regpath = 'input';
for i=1:nargin
    if isa(varargin{i},'datastruct.Dataset')
        data = varargin{i};
    elseif isa(varargin{i},'datastruct.Network')
        net = varargin{i};
    elseif isa(varargin{i},'logical')
        rawZeta = varargin{i};
    % elseif isa(varargin{i},'char')
    %     regpath = varargin{i};
    elseif varargin{i}==floor(varargin{i}) & length(varargin{i}) == 1 & varargin{i} ~= 1 % is integer of length 1 and is not 1
        tmpstraps = varargin{i};
    else
        if isempty(zetavec)
            zetavec = varargin{i};
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
if ~exist('zetavec','var') & strcmpi(regpath,'input')
    estA = clr(data.Y, 'stouffer', 7, 2); %%FOR GENESPIDER DATA
    varargout{1} = estA;
    return
else
    Als = clr(data.Y, 'stouffer', 7, 2); %%FOR GENESPIDER DATA
end

% if strcmpi(regpath,'full')
%     zetavec = abs(Als(:)');
%     zetavec = unique(zetavec);
%     % zetavec = zetavec(zetavec~=0)
% end

zR = []; % zeta range for all bootstraps
Apos = zeros(data.N,data.N,straps);
Alogical = [];
for j = 1:straps
    zetavec = zetavec;
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

if ~rawZeta & strcmpi(regpath,'input')
    zetaRange = [];
    estA = Als;
    zetaRange(1) = min(abs(estA(estA~=0)))-eps;
    zetaRange(2) = max(abs(estA(estA~=0)))+10*eps;
    zR(:,j) = zetaRange;
    varargout{3} = zR;

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

    Als = clr(data.Y, 'stouffer', 7, 2); %%FOR GENESPIDER DATA
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
Asign_frac(isnan(Asign_frac)) = 0;

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