function varargout = RNI(varargin)
% function estA = RNI(data <,net,alpha,regpath>)
%
%   Input Arguments: RNI(data <,net,alpha,zetavec,regpath>)
%   ================
%   data:    datastruct.Dataset
%   net:     datastruct.Network <optional>
%   alpha:   confidense level of inference, can only be a single digit \in [0,1[. default 0.01
%   zetavec: should the confidence values be truncated at some specified values.
%            Will be ignored if regpath = 'full'
%   regpath  {'input','full'}  string to determine if we should try to create a zetavec
%            for the complete regularization path dependant on method, default 'input'.
%            Where 'input' is a zetavec determined by the user
%            and 'full' gives the best estimate of the complete regularization path
%
%   Output Arguments: confA, zetavec
%   =================
%   confA: the confident networks as a 3d array.
%   zetavec: the complete regularization path
%   zetaRange: will return the raw zeta range scale factors
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse input arguments %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawZeta = 0;
regpath = 'input';
alpha = 0.01;
net = [];
for i=1:nargin
    if isa(varargin{i},'datastruct.Dataset')
        data = varargin{i};
    elseif isa(varargin{i},'datastruct.Network')
        net = varargin{i};
    elseif isa(varargin{i},'logical')
        rawZeta = varargin{i};
    elseif isa(varargin{i},'char')
        regpath = varargin{i};
    else
        if length(varargin{i}) > 1
            zetavec = varargin{i};
        else
            alpha = varargin{i};
        end
    end
end

if ~exist('data')
    error('needs a data set')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = data.lambda;
if numel(lambda) == 2
    o = ones(size(data.P));
    [gamma, ps, Nrps] = RInorm(response(data,net)',data.P',diag(lambda(1:length(lambda)/2))*o',diag(lambda(length(lambda)/2+1:end))*o'+eps,alpha);
else
    o = ones(size(data.P),2);
    [gamma, ps, Nrps] = RInorm(response(data,net)',data.P',(diag(lambda(1:length(lambda)/2))'*o)',(diag(lambda(length(lambda)/2+1:end))'*o)'+eps,alpha);
end

if strcmpi(regpath,'full')
    zetavec = abs(gamma(:)');
    zetavec = unique(zetavec);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine how to handle zeta %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~rawZeta & strcmpi(regpath,'input')
    zetaRange = [];
    estA = gamma;
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

for i=1:length(zetavec)
    temp = find(abs(gamma) <= zetavec(i));
    Atmp = gamma;
    Atmp(temp) = 0;

    estA(:,:,i) = Atmp;
end

varargout{1} = estA;

if nargout > 1 & strcmpi(regpath,'input')
    varargout{2} = zetavec;
    varargout{3} = zetaRange;
elseif nargout > 1 & strcmpi(regpath,'full')
    zetavec = (zetavec-zetaRange(1))/delta;
    varargout{2} = zetavec;
    varargout{3} = zetaRange;
end

function varargout = RInorm(varargin)
% RInorm Robust Inference for normally distributed errors
%
%   Input Arguments: RInorm(Phi, Xi, LamPhi, LamXi, alpha)
%   ================
%   Phi    = The regressor matrix, with one variable per column.
%   Xi     = The regressand matrix, with one variable per column.
%   LamPhi = The variance of the regressors (default all 1).
%   LamXi  = The variance of the regressands (default all 1).
%   alpha  = The desired significance level [0,1] (default 0.05).
%
%   Output Arguments: [Gamma, PS, NrPS] 
%   =================
%   Gamma  = The confidence score for practical independence of each 
%            Psi_{ki} = [phi_{1},...,phi_{k-1},phi_{k+1},...,phi_{n},xi_{i}],
%            i.e. gamma_{ik} is the confidence score for practical
%            selectability of regressor phi_{k} with regard to explanation
%            of regressand xi_{i}. The corresponding interaction a_{ik} has
%            the same confidence score for classification as existing.
%            If Xi = [] is given, then the confidence score for practical
%            independence of Phi is returned.
%   PS     = Practically selectable regressors and existing interactions.
%   NrPS   = Number of practically selectable regressors and existing
%            interactions.
%
%   Example
%   =======
%   Assuming the data model Y = -A^{-1}*(P - F) + E, with e_{ij} ~
%   N(0,lamE_{ij}) and f_{ij} ~ N(0,lamF_{ij}), then the confidence scores
%   at significance level 0.01 are obtained by 
%   Gamma = RInorm(Y',-P',LamE,LamF,0.01)
%
%   For further information see 
%   Nordling 2013 Robust inference of gene regulatory networks, PhD thesis,
%   KTH Royal Institute of Technology, Stockholm, Sweden. 
%
%   Information:
%   ============
%   TN's Matlab Toolbox
%   Copyright (c) 2013 Torbjoern Nordling, tn@kth.se
%   All rights reserved.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 0.05; %Default significance level
LamPhi = 1; %Default variance of regressors
LamXi = 1; %Default variance of regressands
Xi = []; %Default is empty regressand matrix
cpus = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1,
    error('NetworkInference:tools:RInorm:InputNumberError','Too few input arguments. At least the regressor matrix must be supplied.');
else
    Phi = varargin{1};
    if nargin > 1;
        Xi = varargin{2};
        if nargin > 2;
            LamPhi = varargin{3};
            if nargin > 3;
                LamXi = varargin{4};
                if nargin > 4;
                    if varargin{5} < 1
                        alpha = varargin{5};
                    else
                        cpus = varargin{5};
                    end
                    if nargin > 5;
                        if varargin{6} < 1
                            alpha = varargin{6};
                        else
                            cpus = varargin{6};
                        end
                    end                    
                end
            end
        end
    end
end

% Checks number of samples
if ~isempty(Xi),
    if size(Phi,1) ~= size(Xi,1),
        error('NetworkInference:tools:RInorm:NumberOfRowsError','The number of row in the regressand matrix does not match the number of rows in the regressor matrix.');
    end
end

% Checks and constructs variance matrices
if any(size(LamPhi) ~= size(Phi)),
    if numel(LamPhi) == 1, %Scalar
        LamPhi = LamPhi.*ones(size(Phi));
    elseif numel(LamPhi) == size(Phi,2), %Vector with as many elements as variables
        LamPhi = ones(size(Phi,1),1)*reshape(LamPhi,1,size(Phi,2));
    else
        error('NetworkInference:tools:RInorm:VarianceMatrixError','The variance matrix of the regressors has the wrong size. It should have the same size as Phi.');
    end
end
if any(LamPhi < 0),
    error('NetworkInference:tools:RInorm:VarianceMatrixError','Some element of the variance matrix of the regressors is negative.');
end
if ~isempty(Xi),
    if any(size(LamXi) ~= size(Xi)),
        if numel(LamXi) == 1, %Scalar
            LamXi = LamXi.*ones(size(Xi));
        elseif numel(LamXi) == size(Xi,2), %Vector with as many elements as variables
            LamXi = ones(size(Xi,1),1)*reshape(LamXi,1,size(Xi,2));
        else
            error('NetworkInference:tools:RInorm:VarianceMatrixError','The variance matrix of the regressands has the wrong size. It should have the same size as Xi.');
        end
    end
    if any(LamXi < 0),
        error('NetworkInference:tools:RInorm:VarianceMatrixError','Some element of the variance matrix of the regressands is negative.');
    end
end

% Checks range of alpha
if alpha > 1 || alpha < 0,
    error('NetworkInference:tools:RInorm:RangeError','The significance level is outside of its range [0,1].');
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE CONFIDENCE SCORE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
Gamma = zeros(size(Phi,2),size(Phi,2));
n = size(Phi,2);
PhiS = Phi./sqrt(chi2inv(1-alpha,prod(size(Phi))).*LamPhi);% Scales all elements of Phi
if ~isempty(Xi), %Practical independence of Psi
    XiS = Xi./sqrt(chi2inv(1-alpha,prod(size(Xi))).*LamXi);% Scales all elements of Xi
    parfor (i=1:size(Xi,2),cpus)
        for j=1:n
            if size(Phi,1) < size(Phi,2),
                Gamma(i,j) = 0;
            else
                Psi = PhiS;
                Psi(:,j) = XiS(:,i);
                SPsi = svd(Psi);
                Gamma(i,j) = SPsi(end);
            end
        end
    end
else %Practical independence of Phi
    if size(Phi,1) < size(Phi,2)
        Gamma = 0;
    else
        Gamma = min(svd(PhiS));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout == 0
    disp(Gamma);
elseif nargout == 1
    varargout{1} = Gamma; 
elseif nargout == 2
    varargout{1} = Gamma; 
    PS = Gamma > 1;
    varargout{2} = PS; 
elseif nargout == 3
    varargout{1} = Gamma; 
    PS = Gamma > 1;
    NrPS = sum(sum(PS));
    varargout{2} = PS; 
    varargout{3} = NrPS; 
end

return
