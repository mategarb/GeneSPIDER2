function A = stabilize(Atilde,varargin)
% Stabilize weights a static network structure
%
%   Input arguments: stabilize(Atilde [,'iaa', iaa, 'sign', sign])
%   ================
%   Atilde:     Network model, with specific desired structure (matrix)
%
%   name : value pairs
%   ================
%   iaa:        wished Interrampatteness, low/high (string) (low)
%               low  => cond < N/2
%               high => cond > N*2
%               N = #nodes
%
%   sign:  if signed structure is to be kept (logical) (false)
%          may not be solveable for some signed structures.
%

%% Parse Input
iaa = 'low';
inputopts = struct('iaa',iaa,'sign',false);
optionNames = fieldnames(inputopts);
nArgs = length(varargin);
if round(nArgs/2) ~= nArgs/2
   error('stabilize needs Name/Value pairs')
end
for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
    inpName = lower(pair{1}); % make case insensitive
    if any(strmatch(inpName,optionNames))
        inputopts.(inpName) = pair{2};
    else
        warning('%s is not a recognized parameter name',inpName)
    end
end

if strcmp(inputopts.iaa,'low')
    Epsilon = -0.01;
    Gamma = -10;
elseif strcmp(inputopts.iaa,'high')
    Epsilon = -0.01;
    Gamma = -100;
    Atilde = 10*Atilde;
end

%% Stabilize matrix
tol = 1e-4;
S = abs(Atilde) < tol;
n = size(Atilde,1);
I = eye(n);
Psi = Gamma*0.9;
Shi = Epsilon*1.1;

if ~inputopts.sign
    try
        cvx_begin sdp quiet
        cvx_precision('low')
        variables D(n,n) g e
        minimize norm(D,'fro')
        subject to
        g*I <= ((Atilde + D) + (Atilde + D)') <= e*I
        Gamma <= g <= Psi
        Shi <= e <= Epsilon
        D(S) == 0
        cvx_end
    catch ME
        fprintf('%s\n',ME.message)
    end
elseif inputopts.sign
    Neg = Atilde < -tol;
    Pos = Atilde > tol;
    try
        cvx_begin sdp
        cvx_precision('low')
        variables D(n,n) g e
        minimize norm(D,'fro')
        subject to
        g*I <= ((Atilde + D) + (Atilde + D)') <= e*I
        Gamma <= g <= Psi
        Shi <= e <= Epsilon
        D(S) == 0
        Atmp = Atilde + D
        Atmp(Pos) >= 0
        Atmp(Neg) <= 0
        cvx_end
    catch ME
        fprintf('%s\n',ME.message)
    end
end

D = full(D);
A = Atilde + full(D);
