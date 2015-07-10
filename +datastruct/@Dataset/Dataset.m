classdef Dataset < datastruct.Exchange
% The Dataset class is used to store a datastruct.Dataset
% complementary to a datastruct.Network
%
%   Input arguments: Dataset([Network ,struct])
%   ================
%   (none) :      Initiate Dataset
%   Network :     A Network object.
%   struct :      A structure with named fields.
%                 If an Experiment is not supplied
%                 then data needs to be added manualy later with the function populate
%
%   ================ functionality
%

    properties (SetAccess = private)
        dataset ='';  % Name of the data set
        network ='';  % Name of complementary network
        P         % True perturbations
        F         % Perturbation noise
        cvP = []; % Covariance of P
        sdP = []; % measurement point variation of P
        Y         % True expression response
        E         % Expression noise
        cvY = []; % Covariance of noisy Y
        sdY = []; % measurement point variation of Y
        lambda    % Noise variance, lambda
        SNR       % Signal to noise ratio, \sigma_N(Y)/\sigma_1(E)
        SNRnu     % Signal to noise ratio, \argmin_i min(norm(\Fi_i)/norm(\nu_i))
        SNRm      % Signal to noise ratio,
        SNRv      % Signal to noise ratio,
        SNRvm     % Signal to noise ratio,
        SNRnum    % Signal to noise ratio,
        etay      % Linear independance from data per sample.
        etau      % Linear independance from data per perturbation.
    end

    properties (SetAccess = public)
        alpha = 0.05; % Confidence
        names = {};
        desc = ''; % description
    end

    properties (Hidden = true)
        N           % # variables in A
        M           % # experiments
        created = struct('creator','','time',now,'id','','nexp','');
        tol = eps;
    end

    methods
        function data = Dataset(varargin)
            warning('off','MATLAB:structOnObject')
            data = data@datastruct.Exchange();
            if nargin > 0
                for i=1:nargin
                    if isa(varargin{i},'datastruct.Network')
                        populate(data,varargin{i});
                    elseif isa(varargin{i},'datastruct.Experiment')
                        experiment = struct(varargin{i});
                        experiment.Y = trueY(varargin{i});
                        populate(data,experiment);
                        data.created.id = num2str(round(cond(data.Y)*10000));
                    elseif isa(varargin{i},'datastruct.Dataset')
                        newdata = struct(varargin{i});
                        populate(data,newdata);
                        data.created.id = num2str(round(cond(data.Y)*10000));
                    elseif isa(varargin{i},'struct')
                        input = varargin{i};
                        populate(data,input);
                        data.created.id = num2str(round(cond(data.Y)*10000));
                    end
                end
                if ispc; data.created.creator = getenv('USERNAME');
                else; data.created.creator = getenv('USER');
                end
                setname(data)
            end
        end

        function M = get.M(data)
            M = size(data.P,2);
        end

        function N = get.N(data)
            N = size(data.P,1);
        end

        function set.lambda(data,lambda)
            if ~isrow(lambda)
                lambda = lambda';
            end
            if prod(size(lambda)) == 1
                lambda = [lambda, 0];
            end
            if prod(size(lambda)) == data.N
                lambda = [lambda, zeros(size(lambda))];
            end

            if ~mod(length(lambda),2) == 0
                error('Something is wrong with the size of lambda. Help!')
            end
            data.lambda = lambda;
        end

        function SNR = get.SNR(data)
            SNR = min(svd(data.Y))/max(svd(data.E));
        end

        function SNR = get.SNRnu(data)
            snr = [];
            for i=1:data.N
                snr(i) = norm(data.Y(i,:))/norm(data.E(i,:));
            end
            SNR = min(snr);
        end

        function SNR = get.SNRnum(data)
            snr = [];
            for i=1:data.N
                snr(i) = norm(data.Y(i,:))/norm(data.E(i,:));
            end
            SNR = mean(snr);
        end

        function SNR = get.SNRm(data)
            alpha = data.alpha;
            sigma = min(svd(data.Y));
            SNR = sigma/sqrt(chi2inv(1-alpha,prod(size(data.P)))*data.lambda(1));
        end

        function SNR = get.SNRv(data)
            alpha = data.alpha;
            for i=1:data.N
                snr(i) = norm(data.Y(i,:))/sqrt(chi2inv(1-alpha,data.M)*data.lambda(1));
            end
            SNR = min(snr);
        end

        function SNR = get.SNRvm(data)
            alpha = data.alpha;
            for i=1:data.N
                snr(i) = norm(data.Y(i,:))/sqrt(chi2inv(1-alpha,data.M)*data.lambda(1));
            end
            SNR = mean(snr);
        end

        function setname(data,varargin)
            if nargin == 1
                namestruct = data.created;
            elseif nargin == 2
                namestruct = varargin{1};
            end
            if ~isa(namestruct,'struct')
                error('input argument must be name/value pairs in struct form')
            end

            namer = data.created;
            optNames = fieldnames(namer);
            inpNames = fieldnames(namestruct);

            for i=1:length(inpNames)
                if any(strcmp(inpNames{i},optNames))
                    namer.(inpNames{i}) = namestruct.(inpNames{i});
                end
            end
            if isempty(data.lambda)
                SNRm = '0';
            else
                SNRm = num2str(round(data.SNRm*1000));
            end
            data.dataset = [namer.creator,'-ID',data.network(regexpi(data.network,'-ID')+3:end),'-D',datestr(namer.time,'yyyymmdd'),'-E',num2str(size(data.P,2)),'-SNR',SNRm,'-IDY',namer.id];
        end

        function names = get.names(data)
            names = data.names;
            if isempty(names)
                for i=1:data.N
                    names{i} = sprintf(['G%0',num2str(floor(log10(data.N))+1),'d'],i);
                end
            else
                names = data.names;
            end
        end

        function varargout = std(data)
            if numel(data.lambda) == 2,
                sdY = data.lambda(1)*ones(size(data.P));
                sdP = data.lambda(2)*ones(size(data.P));
            else
                sdY = data.lambda(1:data.N)'*ones(1,size(data.P,2));
                sdP = data.lambda(data.N+1:end)'*ones(1,size(data.P,2));
            end
            if nargout == 1
               varargout{1} = sdY;
            elseif nargout == 2
               varargout{1} = sdY;
               varargout{2} = sdP;
            end
        end

        function set_new_lambda(data,lambda)
            data.lambda = lambda;
        end

        function scale_lambda_SNRm(data,SNRm)
        % scale the noise variance by setting the theoretic lambda with wished SNRm
            s = min(svd(data.Y));
            lambda = s^2/(chi2inv(1-data.alpha,prod(size(data.P)))*SNRm^2);
            data.lambda = lambda;
            % data.E = sqrt(lambda/preLambda).*data.E;
        end

        function scale_lambda_SNRv(data,SNRv)
        % scale the noise variance by setting the theoretic lambda with wished SNRv
            for i=1:data.N
                lambda(i) = norm(data.Y(i,:))/(chi2inv(1-data.alpha,data.M)*SNRv^2);
            end
            data.lambda = min(lambda);

            % alpha = data.alpha;
            % for i=1:data.N
            %     snr(i) = norm(data.Y(i,:))/sqrt(chi2inv(1-alpha,data.M)*data.lambda(1));
            % end
            % SNR = min(snr);
        end

        function newdata = scaleSNR(data,net,SNR)
        % scales the noise variance to achieve desired SNR. A new dataset is
        % created based on the old data and the new SNR.
        %
        % == Usage ==
        % newdata = scaleSNR(data,net,SNR)
        %
            newdata = datastruct.Dataset(data,net);
            sY = svd(newdata.Y);
            sE = svd(newdata.E);
            scale = 1/SNR*min(sY)/max(sE);
            newdata.lambda = scale^2*newdata.lambda;
            newdata.E = scale*newdata.E;
        end

        function gaussian(data)
        % Generate new gaussian noise matrices E and F with variance lambda for
        % response and/or perturbations.
        %
        % == Usage ==
        % gaussian(data)
        %

            if numel(data.lambda) == 1,
                E = sqrt(data.lambda).*randn(data.N,data.M);
                F = zeros(data.M,data.N);
            elseif numel(data.lambda) == 2,
                E = sqrt(data.lambda(1)).*randn(data.N,data.M);
                F = sqrt(data.lambda(2)).*randn(data.N,data.M);
            elseif numel(data.lambda) == data.N,
                E = sqrt(data.lambda)'*randn(1,data.N);
                F = zeros(data.N,data.M);
            elseif numel(data.lambda) == 2*data.N,
                E = sqrt(data.lambda(1:data.N))'*randn(1,data.N);
                F = sqrt(data.lambda(data.N+1:end))'*randn(1,data.N);
            end
            data.E = E;
            data.F = F;
        end

        function setsdY(data,sdY)
            data.sdY = sdY;
        end

        function setsdP(data,sdP)
            data.sdP = sdP;
        end

        function varargout = cov(data)
        % returns default covariance for the data.
            if numel(data.lambda) == 1
                cvY = data.lambda*eye(data.N);
                cvP = 0*eye(data.N);
            elseif numel(data.lambda) == 2
                cvY = data.lambda(1)*eye(data.N);
                cvP = data.lambda(2)*eye(data.N);
            elseif numel(data.lambda) == data.N
                cvY = diag(data.lambda);
                cvP = 0*eye(size(data.N,1));
            elseif numel(data.lambda) == 2*data.N
                cvY = diag(data.lambda(1:data.N));
                cvP = diag(data.lambda(data.N+1:end));
            end
            if nargout == 1
               varargout{1} = cvY;
            elseif nargout == 2
               varargout{1} = cvY;
               varargout{2} = cvP;
            end
        end

        function setcovY(data,cvY)
           data.cvY = cvY;
        end

        function setcovP(data,cvP)
            data.cvP = cvP;
        end

        function Y = response(data,varargin)
        % Gives the networks response to input from data.

            [n,m] = size(data.P);
            if length(varargin) == 1
                if isa(varargin{1},'datastruct.Network')
                    net = varargin{1};
                    Y = net.G*(data.P-data.F(:,1:m)) + data.E(:,1:m);
                else
                    Y = data.Y+data.E(:,1:m);
                end
            else
                Y = data.Y+data.E(:,1:m);
            end
        end

        function newdata = without(data,i,varargin)
        % creates a Dataset without sample i
        %
        % == Usage ==
        % newdata = without(data,i [,net])
        %           where 'i' is the sample that should be removed
        %           and net is a datastruct.Network.
        %
            M = data.M;

            tmp(1).Y = data.Y;
            tmp(1).P = data.P;
            tmp(1).E = data.E;
            tmp(1).F = data.F;

            tmp(1).Y(:,i) = [];
            tmp(1).P(:,i) = [];
            tmp(1).E(:,i) = [];
            tmp(1).F(:,i) = [];

            sdY = data.sdY;
            if ~isempty(sdY)
                sdY(:,i) = [];
                tmp(1).sdY = sdY;
            end

            sdP = data.sdP;
            if ~isempty(sdP)
                sdP(:,i) = [];
                tmp(1).sdP = sdP;
            end

            if length(varargin) == 1
                newdata = datastruct.Dataset(data,net);
            else
                newdata = datastruct.Dataset(data);
            end

            newdata = populate(newdata,tmp);
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
                varargout{1} = irr;
            end
            if nargout >= 1
                % irr = sum(sum(SIC < 1))/data.N^2;
                varargout{1} = irr;
            end
            if nargout >= 2
                varargout{2} = SIC;
            end
        end

        function eta(data,varargin)
        % calculates eta values for Y and P.
        %
        % eta(data[,net])
        %
            net = [];
            if length(varargin) == 1
                if isa(varargin{1},'datastruct.Network')
                    net = varargin{1};
                end
            end
            Y = response(data,net);
            P = data.P-data.F;

            for i=1:data.M
                Ytemp = Y;
                Ptemp = P;
                Ytemp(:,i) = [];
                Ptemp(:,i) = [];
                ytemp = Y(:,i);
                ptemp = P(:,i);
                etay(i) = sum( abs(Ytemp'*ytemp) );
                etau(i) = sum( abs(Ptemp'*ptemp) );
            end
            data.etay = etay;
            data.etau = etau;
        end

        function w_eta(data,varargin)
        % calculates eta values for Y and P.
        %
        % eta(data[,net])
        %

            net = [];
            if length(varargin) == 1
                if isa(varargin{1},'datastruct.Network')
                    net = varargin{1};
                end
            end
            Y = response(data,net);
            P = data.P-data.F;
            etay = [];
            etau = [];
            for i=1:data.M
                tmpdata = without(data,i);

                [yiU, yiS, yiV] = svd(response(tmpdata,net));
                [uiU, uiS, uiV] = svd(tmpdata.P-tmpdata.F);
                dyiS = diag(yiS);
                duiS = diag(uiS);

                if length(dyiS) < tmpdata.M
                    dyiS = [ dyiS; zeros( (size(tmpdata.N,1)-length(dyiS)),1 ) ];
                    duiS = [ duiS; zeros( (size(tmpdata.N,1)-length(duiS)),1 ) ];
                end

                yv = Y(:, i) / norm( Y(:, i) );
                uv = P(:, i) / norm( P(:, i) );

                yg = dyiS' * abs(yiU' * yv) / sum(dyiS); %yiS(1);
                etay(i) = yg;
                ug = duiS' * abs(uiU' * uv) / sum(duiS); %uiS(1);
                etau(i) = ug;
            end
            data.etay = etay;
            data.etau = etau;
        end

        function included = include(data,varargin)
        % based on eta limit, returns samples ok to include in LOOCO
        %
        % Currently assumes F = 0
        %
        % included = include(data {, etaLim})
        %
        % etaLim is a vector with one element for each limit for Y and P.
        % If one element is given it's used for both Y and P
        %
            if length(varargin) == 1
                etaLim = varargin{1};
            end

            if isempty(data.etay)
                eta(data);
            end

            if ~exist('etaLim','var')
                SY = svd(response(data));
                SP = svd(data.P+data.F);
                included = [and(data.etay >= SY(data.N), data.etau >= SP(data.N))];
            elseif exist('etaLim','var')
                if numel(etaLim) == 1
                    included = [and(data.etay >= etaLim, data.etau >= etaLim)];
                else
                    included = [and(data.etay >= etaLim(1), data.etau >= etaLim(2))];
                end
            end
        end

        function varargout = bootstrap(data,varargin)
        % Bootstraps a new data set from the old data set.
            if length(varargin) == 1
                boots = varargin{1};
            else
                boots = randi([1 data.M],1,data.M);
            end

            tmpdata = datastruct.Dataset();
            tmpdata.populate(data);

            tmp = struct([]);

            tmp(1).P = data.P(:,boots);
            tmp(1).Y = data.Y(:,boots);
            if ~isempty(data.sdY)
                tmp(1).sdY = data.sdY(:,boots);
            end
            if ~isempty(data.sdP)
                tmp(1).sdP = data.sdP(:,boots);
            end
            if ~isempty(data.E)
                tmp(1).E = data.E(:,boots);
            end
            if ~isempty(data.F)
                tmp(1).F = data.F(:,boots);
            end

            tmpdata.populate(tmp);
            varargout{1} = tmpdata;
            if nargout == 2
               varargout{2} = boots;
            end
        end

        function varargout = populate(data,input)
        % populate the Dataset object with matching fields of the input.
        %
        % == Usage ==
        % {data =} populate(data,input)
        %          With input being a struct, datastruct.Dataset,
        %          datastruct.Experiment or datastruct.Network

            if ~isa(input,'struct') && ~isa(input,'datastruct.Dataset') && ~isa(input,'datastruct.Experiment') && ~isa(input,'datastruct.Network')
                error('Needs to be a struct,\n datastruct.Dataset,\n datastruct.Experiment or datastruct.Network class')
            end

            inputnames = fieldnames(input);
            names = fieldnames(data);
            for name = inputnames'
                if any(strmatch(name,names))
                    data.(name{1}) = input.(name{1});
                end
            end

            if nargout > 0
                varargout{1} = data;
            end
        end

        function save(data,varargin)
            save@datastruct.Exchange(data,varargin{:});
        end
    end
    methods (Static)
        function varargout = load(varargin)

            if (length(varargin) == 1) & (exist(varargin{1}) == 7) & (nargout == 0)
                load@datastruct.Exchange(varargin{:});
                return
            elseif (length(varargin) == 1) & (exist(varargin{1}) == 7) & (nargout == 1)
                file_list = load@datastruct.Exchange(varargin{:});
                varargout{1} = file_list;
                return
            end

            [data,dataset] = load@datastruct.Exchange(varargin{:});

            if nargout == 1
                varargout{1} = data;
            elseif nargout == 2
                varargout{1} = data;
                varargout{2} = dataset;
            end
        end
    end
end
