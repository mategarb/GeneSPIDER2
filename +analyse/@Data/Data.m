classdef Data

    properties (SetAccess = private)
        dataset = ''; % identifier for dataset used
        SNR_E         % Signal to noise ratio, \sigma_N(Y)/\sigma_1(E)
        SNR_L         % Signal to noise ratio,
        SNR_e         % Signal to noise ratio, \argmin_i min(norm(\Fi_i)/norm(\nu_i))
        SNR_l         % Signal to noise ratio,
        mean_SNR_e    % Signal to noise ratio,
        mean_SNR_l    % Signal to noise ratio,
        
        % N             % Number of variables in model
        % M             % Number of equations in data

    end

    %  
    
    properties (SetAccess = public)
        tol = eps;         % tolerance level
        alpha = 0.05;      % confidence
    end

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
            analysis.dataset        = analyse.Data.identifier(data);
            analysis.SNR_E          = analyse.Data.calc_SNR_E(data);
            analysis.SNR_L          = analyse.Data.calc_SNR_L(data);
            analysis.SNR_e          = analyse.Data.calc_SNR_e(data);
            analysis.SNR_l          = analyse.Data.calc_SNR_l(data);
            analysis.mean_SNR_e     = analyse.Data.calc_mean_SNR_e(data);
            analysis.mean_SNR_l     = analyse.Data.calc_mean_SNR_l(data);
        end
    end

    methods (Static)
        function dataset = identifier(data)
            if isa(data,'datastruct.Dataset')
                dataset = data.dataset;
            else
                dataset = '';
            end
        end

        % function [N,M] = size(data)
        %     N = data.N;
        %     M = data.M;
        % end

        function SNR = calc_SNR_E(data)
            SNR = min(svd(data.Y))/max(svd(data.E));
        end

        function SNR = calc_SNR_e(data)
            snr = [];
            for i=1:data.N
                snr(i) = norm(data.Y(i,:))/norm(data.E(i,:));
            end
            SNR = min(snr);
        end

        function SNR = calc_mean_SNR_e(data)
            snr = [];
            for i=1:data.N
                snr(i) = norm(data.Y(i,:))/norm(data.E(i,:));
            end
            SNR = mean(snr);
        end

        function SNR = calc_SNR_L(data)
            alpha = data.alpha;
            sigma = min(svd(data.Y));
            SNR = sigma/sqrt(chi2inv(1-alpha,prod(size(data.P)))*data.lambda(1));
        end

        function SNR = calc_SNR_l(data)
            alpha = data.alpha;
            for i=1:data.N
                snr(i) = norm(data.Y(i,:))/sqrt(chi2inv(1-alpha,data.M)*data.lambda(1));
            end
            SNR = min(snr);
        end

        function SNR = calc_mean_SNR_l(data)
            alpha = data.alpha;
            for i=1:data.N
                snr(i) = norm(data.Y(i,:))/sqrt(chi2inv(1-alpha,data.M)*data.lambda(1));
            end
            SNR = mean(snr);
        end
        
        function varargout = irrepresentability(data,net)
        % Calculates the irrepresentability of the data set for inference with
        % LASSO
            Y = data.Y;
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
