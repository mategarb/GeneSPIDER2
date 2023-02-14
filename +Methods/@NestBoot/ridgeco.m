function varargout = ridgeco(varargin)
%
% function estA = Boridge(data,net,zetavec[,alpha,rawZeta,straps])
%
%   Input Arguments: Boridge(data,net,zetavec[,alpha,rawZeta,straps])
%   ================
%   data:    datastruct.Dataset
%   net:     datastruct.Network
%   zetavec: method parameter for tuning the network fitting.
%   alpha:   The elasticnet mixing parameter, with 0 < alpha <= 1. (default = 1, Lasso)
%            Currently alpha < 0.01 is not reliable, unless you
%            supply your own zeta sequence. zetavec needs to be set first
%   rawZeta: logical to determine if the zeta values should be
%            converted.  default = false
%   straps:  integer, number of bootstrap runs. must be > 1, default = 100
%
%   Output Arguments: estA
%   =================
%   estA: the estimated networks as a 3d array.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse input arguments %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawZeta = 0;
tmpzetas = [];
net = [];
tmpstraps = 1;
regpath = 'input';
for i = 1:nargin
    if isa(varargin{i},'datastruct.Dataset')
        data = varargin{i};
    elseif isa(varargin{i},'datastruct.Network')
        net = varargin{i};
    elseif isa(varargin{i},'logical')
        rawZeta = varargin{i};
    elseif all(varargin{i} == floor(varargin{i})) && all(length(varargin{i}) == 1) && all(varargin{i} ~= 1) % is integer of length 1 and is not 1
        tmpstraps = varargin{i};
    else
        if isempty(tmpzetas)
            tmpzetas = varargin{i};
        %else
        %    alpha = varargin{i};
        end
    end
end

if ~exist('data','var')
    error('needs a data set')
end

if tmpstraps == 1
    straps = 100;
else
    straps = tmpstraps;
end

zR = zeros(2,straps); % zeta range for all bootstraps
nlinksBo = zeros(length(tmpzetas),straps);
%Apos = zeros(data.N,data.N,straps);
Alogical = [];
for j = 1:straps
    zetavec = tmpzetas;
    bdata = bootstrap(data);

    reps = 1;
    while rank(bdata.P) < min(bdata.N,bdata.M)
        bdata = bootstrap(data);
        if ~mod(reps,10000)
            disp('number of boot tries')
            disp(reps)
        end
        reps = reps + 1;
    end

    %% Run
    Afit = zeros(size(data.P,1),size(data.P,1),length(zetavec));
    for i = 1:size(bdata.P,1)
        %[b, stats] = lasso(response(bdata,net)', -bdata.P(i,:)', 'Lambda', zetavec,'Alpha',eps(1));
        b = ridge(-data.P(i,:)',response(data,net)', zetavec, 1);
        Afit(i,:,:) = b(:,:);
    end
    %for i = 1:size(bdata.P,1)
    %    fit = glmnet(response(bdata,net)',-bdata.P(i,:)','gaussian',glmnetSet(struct('lambda',zetavec,'alpha',alpha)));
    %    try
    %        Afit(i,:,:) = fit.beta(:,:);
    %    catch
    %        save(bdata,'../debug/')
    %        return
    %    end
    %end
    %Afit(:, :, :) = Afit(:, :, end:-1:1); % Glmnet reverses the order. Need to undo.
    
    TSS = sum((bdata.Y-mean(bdata.Y)).^2);

    for i = 1:size(Afit,3)
        yfit = -pinv(Afit(:,:,i))*bdata.P;
        RSS = sum((bdata.Y-yfit).^2);
        Rsquared(i) = 1 - RSS/TSS;
    end

    % select the highest scoring dataset 
    [~,select] = max(Rsquared);
    L2_pen = flip(zetavec);
    L2_pen = L2_pen(select);
    % and set this to the new Als
    Als = Afit(:,:,select);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine how to handle zeta for cutoff %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    zetavec=zetavec;
    zetaRange(2) = max(zetavec);
    zetaRange(1) = min(zetavec);
    delta = zetaRange(2)-zetaRange(1);
elseif rawZeta & strcmpi(regpath,'full')
    zetavec=zetavec;
    zetaRange(2) = 1;
    zetaRange(1) = 0;
    delta = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop over each sparsity and each matrix from the fit 
% and remove values smaller than zeta
for i=1:length(zetavec)
    temp = find(abs(Als) <= zetavec(i));
    Atmp = Als;
    Atmp(temp) = 0;

    estA(:,:,i) = Atmp;
end

    if isempty(Alogical)
        Alogical = double(logical(Afit));
        Aposs = double(sign(Afit));
        Aposs(Aposs < 0) = 0;
        Aneg = double(sign(Afit));
        Aneg(Aneg > 0) = 0;
        Aneg = abs(Aneg);
        nlinksBo(:,j) = squeeze(sum(sum(double(logical(Afit)))));
    else
        Alogical = Alogical + double(logical(Afit));
        tmp = double(sign(Afit));
        tmp(tmp < 0) = 0;
        Aposs = Aposs + double(tmp);
        tmp = double(sign(Afit));
        tmp(tmp > 0) = 0;
        Aneg = Aneg + double(abs(tmp));
        nlinksBo(:,j) = squeeze(sum(sum(double(logical(Afit)))));
    end
end

Alogical = Alogical/straps;
Afrac = Alogical;
Alogical(Alogical < 1) = 0;

Asign_Max = zeros(data.N,data.N,size(Aposs,3));
for i = 1:size(Aposs,3)
    tmp = cat(3,Aposs,Aneg);
    Asign_Max(:,:,i) = max(tmp,[],3);
end
Asign_Max = (Asign_Max.*Afrac)/straps;

Aposs_frac = Aposs./(straps*Afrac);
Asign_frac = 2*Aposs_frac-1;
%Asign_frac(isnan(Asign_frac)) = 0;

if nargout > 0
    varargout{1} = Alogical;
end

if nargout > 1 % Link support [0,1]
    varargout{2} = Afrac;
end

if nargout > 2 % Agnostic sign support
    varargout{3} = Asign_Max;
end

if nargout > 3 % Explicit sign support -1 for 100% negative
    varargout{4} = Asign_frac;
end

if nargout > 4 % # positive links
    varargout{5} = cat(3,Aposs,Aneg);
end

if nargout > 5
    varargout{6} = nlinksBo;
end

if nargout > 6
    varargout{7} = zetavec;
end

return
