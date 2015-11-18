function varargout = Bolasso(varargin)
% A function wrapper that implements bootstrap lasso, the lasso algorithm used is the glmnet algorithm:
% (Glmnet for Matlab (2013) Qian, J., Hastie, T., Friedman, J., Tibshirani, R. and Simon, N.
% http://www.stanford.edu/~hastie/glmnet_matlab/)
%
% function estA = Bolasso(data,net,zetavec[,alpha,rawZeta,straps])
%
%   Input Arguments: Bolasso(data,net,zetavec[,alpha,rawZeta,straps])
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

zR = []; % zeta range for all bootstraps
Apos = zeros(data.N,data.N,straps);
Alogical = [];
for j=1:straps
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

    %% Determine how to handle zeta %%

    if ~rawZeta
        zetaRange = [];
        tol = 1e-6;
        zmax = 1;
        zmin = 0;
        % find zero network
        estA = Methods.Glmnet(bdata,net,zmax,logical(1));
        while nnz(estA) > 0
            tmp = zmax;
            zmax = zmax*2;
            estA = Methods.Glmnet(bdata,net,zmax,logical(1));
        end
        % refine
        while zmax-zmin > tol
            i = (zmax + zmin) * 0.5;
            estA = Methods.Glmnet(bdata,net,i,logical(1));
            if nnz(estA) == 0
                zmax = i;
            else
                zmin = i;
            end
        end

        zetaRange(1) = 0;
        zetaRange(2) = zmax;
        zR(:,j) = zetaRange;
        varargout{3} = zR;
        % Convert to interval.
        delta = zetaRange(2)-zetaRange(1);
        zetavec = zetavec*delta + zetaRange(1);
    end

    %% Run

    for i = 1:size(bdata.P,1)
        fit = glmnet(response(bdata,net)',-bdata.P(i,:)','gaussian',glmnetSet(struct('lambda',zetavec,'alpha',alpha)));
        Afit(i,:,:) = fit.beta(:,:);
    end
    Afit(:, :, :) = Afit(:, :, end:-1:1); % Glmnet reverses the order. Need to undo.
    if isempty(Alogical)
        Alogical = double(logical(Afit));
        Aposs = double(sign(Afit));
        Aposs(Aposs < 0) = 0;
        Aneg = double(sign(Afit));
        Aneg(Aneg > 0) = 0;
        Aneg=abs(Aneg);
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

for i=1:size(Aposs,3)
    tmp=cat(3,Aposs,Aneg);
    Asign_Max(:,:,i)=max(tmp,[],3);
end
Asign_Max=(Asign_Max.*Afrac)/straps;

Aposs_frac = Aposs./(straps*Afrac);
Asign_frac = 2*Aposs_frac-1;

if nargout > 0
    varargout{1} = nlinksBo;
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
