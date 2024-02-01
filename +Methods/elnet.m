function varargout = elnet(data, zetavec, rawZeta, alpha)
% 
% This is a wrapper function for elastic net (lasso with alpha 0.7) from matlab toolbox
%
% function estA = elasticNet(data,net,zetavec,rawZeta)
%
%   Input Arguments: elasticNet(data,net,zetavec,rawZeta)
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
%   fit: The glmnet algorithm output data, outputs the last fit if
%        a range of sparsity penalties are supplied.
%   zetaRange: the sparsity range used when a scaled zeta are used (default).
%

arguments
    data (1,1) {mustBeA(data, 'datastruct.Dataset')}
    zetavec (1,:) {mustBeNumeric}
    rawZeta (1,1) {mustBeNumericOrLogical} = 0;
    alpha (1,1) {mustBeNumeric,mustBeReal} = 0.7; % default 0.7, can be tunable
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse input arguments  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('data','var')
    error('Data set is missing')
end

%% Determine how to handle zeta %%

if ~rawZeta
    zetaRange = [];
    tol = 1e-6;
    zmax = 1;
    zmin = 0;
    % find zero network
    estA = Methods.elnet(data,zmax,true);
    while nnz(estA) > 0
        zmax = zmax*2;
        estA = Methods.elnet(data,zmax,true);
    end
    % refine
    while zmax-zmin > tol
        i = (zmax + zmin) * 0.5;
        estA = Methods.elnet(data,i,true);
        if nnz(estA) == 0
            zmax = i;
        else
            zmin = i;
        end
    end

    zetaRange(1) = 0;
    zetaRange(2) = zmax;
    varargout{3} = zetaRange;
    % Convert to interval
    delta = zetaRange(2) - zetaRange(1);
    zetavec = zetavec*delta + zetaRange(1);
end

%% Run
Afit = zeros(size(data.P,1),size(data.P,1),length(zetavec));
for i = 1:size(data.P,1)

   [b,fitinfo] = lasso(data.Y', -data.P(i,:)', 'Lambda', zetavec, 'Alpha', alpha);
   Afit(i,:,:) = b(:,:);

end

    varargout{1} = Afit;

if ~rawZeta
    varargout{2} = zetaRange;
    varargout{3} = fitinfo;
else
    varargout{2} = fitinfo;
end
