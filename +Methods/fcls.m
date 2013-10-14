function varargout = fcls(varargin)
% function Aest = fcls(data,Ainit)
%
% Solve CLS with an unconstrained Newton method. Really it's an
% equality-constrained problem, but I've eliminated the constraints with an
% affine function.
%
%   Input Arguments: mle(data,initA)
%   ================
%   data   = The data in a struct array. Must contain the
%            covariance of 'nY' and 'P' as 'cvY' and 'cvP'
%   Ainit  = Precomputed network matrices as 3D array.
%
%   Output Arguments: Aest
%   =================
%   Aest = The estimated network matrices as cell array.

    net = [];
    verbose = false;
    for i=1:nargin
        if isa(varargin{i},'GeneSpider.Dataset')
            data = varargin{i};
        elseif isa(varargin{i},'GeneSpider.Network')
            net = varargin{i};
        elseif isa(varargin{i},'double')
            Ainit = varargin{i};
        elseif isa(varargin{i},'logical')
            verbose = varargin{i};
        end
    end

    if ~exist('data')
        error('needs a data set')
    end

    %% Parse data to use
    P = data.P;
    Y = response(data,net);

    [cvY,cvP] = cov(data);

    if ~isempty(data.cvY)
        cvY = data.cvY;
    end

    if ~isempty(data.cvP)
        cvP = data.cvP;
    end
    
    e = 1e-6;

    for j=1:size(Ainit,3)
        A = Ainit(:,:,j);
        [H, D, H] = svd(A * data.cvY * A' + data.cvP + e*eye(size(A)));
        dD = diag(D);
        R = H * sqrt(diag( 1 ./dD ));
        results = tools.fcls(A,response(data,net),data.P,R,verbose);

        Aest(:,:,j) = results.A;
    end

    varargout{1} = Aest;
end