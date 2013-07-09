classdef Experiment  < hgsetget
% Experiment by doing perturbations on network A
%
%  This class will automate the process of simulating pertubration
%  experiments in linear systems
%  Data model: Y = -pinv(A)*(P-F) + E, i.e. Y = trueY + E, P = Ptrue + F
%
%  possible settings: set(experiment,'<name>',value)
%
%  'lambda'           Variance of the Gaussian noise that is added ([1 0]).
%  'alpha'            Significance level, 1>alpha>0 (0.05)
%  'nExp'             Number of experiments that can be performed, >0 (unlimited)
%  'mag'              Initial perturbation magnitude
%  'SignalThreshold'  >0 (1)
%  'tol'              tolerance used for calculations (1e-6)
%  'nnzP'             number of non-zero elements for each perturbation sample (size(A,1))
%

    properties (SetAccess = protected)
        A      % True network model
        G      % True static gain matrix
        P = [] % Designed perturbations
        Yinit = []; % initial observed responses
        E = [] % Expression Noise
        F = [] % Perturbation noise
    end

    properties (SetAccess = public)
        lambda = [1 0];       % Variance of the Gaussian noise that is added (1)
        alpha  = 0.05;        % Significance level, 1>alpha>0 (0.05)
        nExp   = inf;         % Number of experiments that can be performed, >0 (unlimited)
        mag    = 1;           % Initial perturbation magnitude
        SignalThreshold = 1;  % >0 (1)
        tol = eps;            % tolerance used for calculations
        nnzP                  % number of non-zero elements for each perturbation sample
                              % design = 'BCSE';      % perturbation design algorithm
        SNR = [];             % Signal to noise ratio
    end

    properties (Hidden = true, SetAccess = protected)
        M             % # variables in A
    end

    methods
        function data = Experiment(varargin)
        % Initiate data to experiment on.
        %
        %   Input Arguments: Experiment(A[, P, Y])
        %   ================
        %   (none) :    Initiate object for later use (a warning will be issued,
        %                                              to supply a network before using object).
        %   A :         The 'true' weighted adjacency matrix of the network.
        %   P :         Previously derived perturbations
        %   Y :         Previously observed responses
        %

            if nargin == 0
                warning('Make sure to supply a Network before further experimentation!')
            end
            if nargin > 0
                data.A = varargin{1};
                % data.lambda(min(min(data.A)));
            end
            if nargin > 1
                if isa(varargin{2},'double')
                    data.P = varargin{2};
                else
                    error('Perturbatin needs to be of type double')
                end
            end

            if nargin > 2
                if isa(varargin{3},'double')
                    data.Yinit = varargin{3};
                else
                    error('initial response needs to be of type double')
                end
            end
        end

        function set.A(data,net)
            if isempty(data.A)
                if isa(net,'double')
                    data.A = net;
                elseif isa(net,'Network')
                    data.A = net.A;
                else
                    error('Network input argument must be a of double array or Network type')
                end
                data.G = -pinv(data.A);
                data.P = [1; zeros(size(data.A,1)-1,1)];
                if isempty(data.nnzP)
                    data.nnzP = data.M;
                end
            else
                warning('True A already set, will not change it!')
            end
        end

        function M = get.M(data)
            M = min(size(data.A));
        end

        function set.alpha(data,SignificanceLevel)
            if SignificanceLevel > 1 || SignificanceLevel < 0
                error('SignificanceLevel must be in the span (0,1)')
            end
            data.alpha = SignificanceLevel;
        end

        function Pinit(data,P)
        % Will set an inital perturbation, overwrites current perturbation
            data.P = P;
        end

        function set.Yinit(data,Yinit)
            data.Yinit = Yinit;
        end

        function set.lambda(data,lambda)
            if ~isrow(lambda)
                lambda = lambda';
            end
            if prod(size(lambda)) == 1
                lambda = [lambda, 0];
            end
            if prod(size(lambda)) == size(data.P,1)
                lambda = [lambda, zeros(size(lambda))];
            end

            if ~mod(length(lambda),2) == 0
                error('Something is wrong with the size of lambda. Help!')
            end
            data.lambda = lambda;
        end

        function set.SignalThreshold(data,SignalThreshold)
            if SignalThreshold <= 0
                error('SignalThreshold must be > 0')
            end
            data.SignalThreshold = SignalThreshold;
        end

        function set.nExp(data,N)
            if rem(N,1) ~= 0 || N < 1
                error('# experiments must be a positive integer')
            end
            data.nExp = N;
        end

        function TF = terminate(data,varargin)
        % Checks whether the stop conditions specified are met.
        %
        %   Input Arguments: terminating(data [,condition])
        %   ================
        %   data       : Experiment object
        %   condition  : use stop conditions
        %                (none) = use only # of samples. This is always on, default to inf.
        %                'ST' = SignalThreshold
            TF = false;
            if size(data.P,2) >= data.nExp
                TF = true;
            end
            if nargin == 2;
                condition = varargin{1};
                if strcmp(condition,'ST') % Signal Threshold
                    s = svd(noiseY(data));
                    if (min(s) > data.SignalThreshold) && (size(data.P,2) >= data.M)
                        TF = true;
                    end
                elseif strcmp(condition,'SCT') % Scaled Threshold
                    k = size(data.P,2);
                    s = svd(noiseY(data)./sqrt(chi2inv(1-data.alpha,data.M*k).*var(data)));
                    if (min(s) > data.SignalThreshold) && (size(data.P,2) >= data.M)
                        TF = true;
                    end
                end
            end
        end

        function Y = trueY(data)
        % generates expression without noise for complete data set
        % Y = trueY(data)
            pre = size(data.Yinit,2) + 1;
            Y = [data.Yinit, data.G*data.P(:,pre:end)];
        end

        function Y = noiseY(data)
        % generates expression with noise for complete data set
        %
        % Y = noiseY(data)
        % also calculate noise if not enough of if has been generated
        %
            pre = size(data.Yinit,2) + 1;
            while size(data.E,2) < size(data.P,2)
                gaussian(data);
            end
            Y = [data.Yinit, data.G*(data.P(:,pre:end)-data.F(:,pre:size(data.P,2))) + data.E(:,pre:size(data.P,2))];
        end

        function P = noiseP(data)
        % Generates expression with noise for complete data set
        %
        % P = noiseP(data)
        % also calculate noise if not enough of if has been generated
            while size(data.F,2) < size(data.P,2)
                gaussian(data);
            end
            P = data.P-data.F(:,1:size(data.P,2));
        end

        function vY = var(data)
        % claculate data point-vise variance of Y
            if numel(data.lambda) == 2,
                vY = data.lambda(1)*ones(size(data.P));
            else
                vY = data.lambda(1:data.M)'*ones(1,size(data.P,2));
            end
        end

        function sparse(data)
        % make each purturbation sufficiently sparse naively
            p = data.P(:,end);
            nZero = data.M-data.nnzP;
            [sortedValues,sortIndex] = sort(abs(p(:)));
            minIndex = sortIndex(1:nZero);
            p(minIndex) = 0;
            data.P(:,end) = p;
        end

        function gaussian(data)
        % generate gaussian noise with variance lambda for response and/or
        % perturbations.
        %
            if numel(data.lambda) == 1,
                data.E = [data.E sqrt(data.lambda).*randn(data.M,1)];
                data.F = [data.F zeros(data.M,1)];
            elseif numel(data.lambda) == 2,
                data.E = [data.E sqrt(data.lambda(1)).*randn(data.M,1)];
                data.F = [data.F sqrt(data.lambda(2)).*randn(data.M,1)];
            elseif numel(data.lambda) == data.M,
                data.E = [data.E sqrt(data.lambda)'.*randn(data.M,1)];
                data.F = [data.F zeros(data.M,1)];
            elseif numel(data.lambda) == 2*data.M,
                data.E = [data.E sqrt(data.lambda(1:data.M))'.*randn(data.M,1)];
                data.F = [data.F sqrt(data.lambda(data.M+1:end))'.*randn(data.M,1)];
            end
        end

        function SNR = get.SNR(data)
            if isempty(data.SNR)
                SNR = min(svd(trueY(data)))/max(svd(data.E));
            else
                SNR = data.SNR;
            end
        end

        function scaleSNR(data,SNR)
        % scales the noise variance to achieve desired SNR, this function sets
        % lambda and change E.
        %
        % == Usage ==
        % scaleSNR(data,SNR)
        %
            noiseY(data);
            sY = svd(trueY(data));
            sE = svd(data.E);
            scale = 1/SNR*min(sY)/max(sE);
            data.lambda = scale^2*data.lambda;
            data.E = scale*data.E;
        end

        function scaleSNRm(data,SNR)
        % scales the noise variance to achieve desired SNR, this function sets
        % lambda and changes E.
        %
        % == Usage ==
        % scaleSNR(data,SNR)
        %
            sY = svd(trueY(data));
            data.lambda = min(sY)^2/(chi2inv(data.alpha,prod(size(data.P)))*SNR^2);
            data.E = [];
            noiseY(data);
        end

        function SVDE(data)
        % SVD     - Linear combination of all directions scaled by their SVs plus a
        %           new orthogonal dimension while needed. Currently assumes F = 0
        %           <=> Ptrue = P.
            k = size(data.P,2);
            if k+1 <= data.M
                newdir = GramSchmidtOrth(data.P(:,1:k),k+1);
                data.P(:,k+1) = data.mag*newdir(:,k+1);
            end
            if k > 0,
                [u,s,v] = svd(noiseY(data));
                if min(diag(s)) < data.SignalThreshold,
                    for j=1:min(size(s)),
                        if s(j,j) > tol,
                            data.P(:,k+1) = data.P(:,k+1) + data.SignalThreshold/s(j,j)*data.P(:,1:k)*v(:,j);
                        end
                    end
                end
            end
        end

        function BCSE(data)
        % BelowCS - Linear combination of all directions scaled by their SVs that
        %           are below the threshold, unless P is singular and a new orthogonal
        %           dimesion is perturbed. Each element of Ytilde is scalled
        %           Ytilde_ij = Y_ij / sqrt(chi^-2(alpha,n*m) var(Y_ij)) such that
        %           the smallest singular value is equal to the confidence score
        %           defined in Nordling 2013 PhD thesis ch. 7.3. var(Y_ij) is the
        %           variance of each element of Y, which is assumed independent and
        %           normally distributed with this variance. Note lambda vector is
        %           used to generate noise of new experiment, while Lambda_ij is
        %           used for old ones. Ysvd and Psvd returns the singular values of
        %           the scaled variables Ytilde and Ptilde.

            k = size(data.P,2);
            if k+1 <= data.M
                newdir = GramSchmidtOrth(data.P(:,1:k),k+1);
                data.P(:,k+1) = data.mag*newdir(:,k+1);
            else
                [uu,su,vu] = svd(noiseY(data));
                [u,s,v] = svd(noiseY(data)./sqrt(chi2inv(1-data.alpha,data.M*k).*var(data)));
                data.P(:,k+1) = zeros(data.M,1);
                for j = 1:data.M,
                    if s(j,j) < data.SignalThreshold,
                        data.P(:,k+1) = data.P(:,k+1) + data.SignalThreshold/s(j,j)*data.P(:,1:k)*vu(:,j);
                    end
                end
                if ~any(diag(s) < data.SignalThreshold)
                    data.P(:,k+1) = data.SignalThreshold/s(j,j)*data.P(:,1:k)*vu(:,j);
                end

            end
        end

        function BSVE(data)
        % BelowSV - Linear combination of all directions scaled by their SVs that
        %           are below the threshold, unless P is singular and a new orthogonal
        %           dimesion is perturbed. Currently assumes F = 0 <=> Ptrue = P.
            k = size(data.P,2);
            if k+1 <= data.M
                newdir = GramSchmidtOrth(data.P(:,1:k),k+1);
                data.P(:,k+1) = data.mag*newdir(:,k+1);
            else
                [u,s,v] = svd(noiseY(data));
                for j = 1:data.M,
                    if s(j,j) < data.SignalThreshold,
                        data.P(:,k+1) = data.SignalThreshold/s(j,j)*data.P(:,1:k)*v(:,j);
                    end
                end
            end
        end

        function PSVE(data)
        % PureSV  - Perturbation in direction corresponding to smallest SV of Y,
        %           unless P is singular and a new orthogonal dimesion is
        %           perturbed. Currently assumes F = 0 <=> Ptrue = P.
            k = size(data.P,2);
            if k+1 <= data.M
                newdir = GramSchmidtOrth(data.P(:,1:k),k+1);
                data.P(:,k+1) = data.mag*newdir(:,k+1);
            else
                [u,s,v] = svd(noiseY(data));
                if min(diag(s)) < data.SignalThreshold,
                    data.P(:,k+1) = data.SignalThreshold/s(data.M,data.M)*data.P(:,1:k)*v(:,data.M);
                end
            end
        end

        function RANDPE(data)
        % RandomP - Random perturbations until stop condition is reached, used in PNAS2011
        %           The perturbation is scaled such that the energy is equal to
        %           the same perturbation of another set of experiments that have been
        %           design previously and is nonzero.
            temp = randn(size(data.G,2),1);
            if k+1 <= data.M,
                temp = temp + GramSchmidtOrth(data.P(:,1:k),k+1);
            end

            if any(abs(data.P(:,k+1)) > tol ),
                data.P(:,k+1) = temp.*(norm(Pall(:,k+1,i-1))/norm(temp));
            else
                data.P(:,k+1) = data.mag*temp./norm(temp);
            end
        end

    end
end
