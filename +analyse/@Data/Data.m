classdef Data < analyse.DataModel
% Data calculates data properties for the supplied data set.
%
%
%
    properties (SetAccess = private)
        dataset = '';     % identifier for dataset used
        SNR_Phi_true      % Signal to noise ratio, \sigma_N(Y)/\sigma_1(E)
        SNR_Phi_gauss     % Signal to noise ratio,
        SNR_phi_true      % Signal to noise ratio, \argmin_i min(norm(\Fi_i)/norm(\nu_i))
        SNR_phi_gauss     % Signal to noise ratio,
    end

    %

    methods
        function analysis = Data(data,varargin)

            if nargin > 1
                ind = find(strcmpi('tol',varargin));
                if ind
                    analysis.tol = varargin{ind+1};
                    varargin(ind,ind+1) = [];
                end
            end

            analysis = analyse_data(analysis,data,varargin{:});

        end

        function analysis = analyse_data(analysis,data,varargin);
            analysis.dataset       = analyse.Data.identifier(data);
            analysis.SNR_Phi_true  = analyse.Data.calc_SNR_Phi_true(data);
            analysis.SNR_Phi_gauss = analyse.Data.calc_SNR_Phi_gauss(data);
            analysis.SNR_phi_true  = min(analyse.Data.calc_SNR_phi_true(data));
            analysis.SNR_phi_gauss = min(analyse.Data.calc_SNR_phi_gauss(data));
        end
    end

    methods (Static)
        function val = alpha(varargin)
        % significance level (default = 0.01)
            val = alpha@analyse.DataModel(varargin{:});
        end

        function val = type(varargin)
        % 'directed' or 'undirected'
            val = type@analyse.DataModel(varargin{:});
        end

        function val = tol(varargin)
        % if a tolernace value for computations is needed it can be set here.
            val = tol@analyse.DataModel(varargin{:});
        end

        function dataset = identifier(data)
            if isa(data,'datastruct.Dataset')
                dataset = data.dataset;
            else
                dataset = '';
            end
        end

        function SNR = calc_SNR_Phi_true(data)
            SNR = min(svd(true_response(data)))/max(svd(data.E));
        end

        function snr = calc_SNR_phi_true(data)
            snr = [];
            X = true_response(data);
            for i=1:data.N
                snr(i) = norm(X(i,:))/norm(data.E(i,:));
            end
        end

        function SNR = calc_SNR_Phi_gauss(data)
            alpha = analyse.Data.alpha;
            sigma = min(svd(response(data)));
            SNR = sigma/sqrt(chi2inv(1-alpha,prod(size(data.P)))*data.lambda(1));
        end

        function SNR = calc_SNR_phi_gauss(data)
            alpha = analyse.Data.alpha;
            Y = response(data);
            for i=1:data.N
                SNR(i) = norm(Y(i,:))/sqrt(chi2inv(1-alpha,data.M)*data.lambda(1));
            end

        end

        function lambda = scale_lambda_SNR_L(data,SNR_L)
        % scale the noise variance by setting the theoretic lambda with wished SNR_L
            alpha = analyse.Data.alpha;
            s = min(svd(true_response(data)));
            lambda = s^2/(chi2inv(1-alpha,prod(size(data.P)))*SNR_L^2);
            % data.E = sqrt(lambda/preLambda).*data.E;
        end

        function lambda = scale_lambda_SNRv(data,SNRv)
        % scale the noise variance by setting the theoretic lambda with wished SNRv
            alpha = analyse.Data.alpha;
            Y = response(data);
            for i=1:data.N
                lambda(i) = norm(Y(i,:))/(chi2inv(1-alpha,data.M)*SNRv^2);
            end
        end

        function varargout = irrepresentability(data,net)
        % Calculates the irrepresentable condition of the data set for inference with LASSO.
            Y = response(data,net);
            Phi = Y';
            for i = 1:net.N
                Phiz = Phi(:,net.A(i,:) == 0);
                Phipc = Phi(:,net.A(i,:) ~= 0);

                sic = abs(Phiz'*Phipc*pinv(Phipc'*Phipc)*sign(net.A(i,net.A(i,:)~=0))');
                k = 1;
                for j=1:net.N
                    if net.A(i,j) == 0
                        SIC(i,j) = sic(k);
                        k = k + 1;
                    end
                end
            end
            % SIC = SIC;
            irr = min(1-SIC(~logical(net)));
            if nargout == 0
                % irr = sum(sum(SIC < 1))/data.N^2;
                disp(irr)
            end
            if nargout >= 1
                % irr = sum(sum(SIC < 1))/data.N^2;
                varargout{1} = irr;
            end
            if nargout >= 2
                varargout{2} = SIC;
            end
        end

    end
end
