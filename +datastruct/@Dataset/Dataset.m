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
        P         % Observed/assumed perturbations
        F         % Perturbation noise
        cvP = []; % Covariance of P
        sdP = []; % measurement point variation of P
        Y         % Observed expression response
        E         % Expression noise
        cvY = []; % Covariance of noisy Y
        sdY = []; % measurement point variation of Y
        lambda    % Noise variance, lambda
        SNR_L     % Signal to noise ratio,
    end

    properties (SetAccess = public)
        names = {};
        description = '';
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
                        data.created.id = num2str(round(cond(data.Y-data.E)*10000));
                    elseif isa(varargin{i},'datastruct.Dataset')
                        newdata = struct(varargin{i});
                        populate(data,newdata);
                        data.created.id = num2str(round(cond(data.Y-data.E)*10000));
                    elseif isa(varargin{i},'struct')
                        input = varargin{i};
                        populate(data,input);
                        data.created.id = num2str(round(cond(data.Y-data.E)*10000));
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

        function SNR = get.SNR_L(data)
            alpha = 0.01;
            sigma = min(svd(true_response(data)));
            SNR = sigma/sqrt(chi2inv(1-alpha,prod(size(data.P)))*data.lambda(1));
        end

        function p = Phi(data)
        % Y'
            p = data.Y';
        end

        function x = Xi(data)
        % -P'
            x = -data.P';
        end

        function u = Upsilon(data)
        % E'
            u = data.E';
        end

        function p = Pi(data)
        % -F'
            p = -data.F';
        end

        function p = Psi(data)
        % [Phi Xi]
            p = [data.Phi data.Xi];
        end

        function o = Omicron(data)
        % [Upsilon Pi]
            o = [data.Upsilon, data.Pi];
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
                SNR_L = '0';
            else
                SNR_L = num2str(round(data.SNR_L*1000));
            end
            data.dataset = [namer.creator,'-ID',data.network(regexpi(data.network,'-ID')+3:end),'-D',datestr(namer.time,'yyyymmdd'),'-N',num2str(size(data.P,1)),'-E',num2str(size(data.P,2)),'-SNR',SNR_L,'-IDY',namer.id];
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
                sdY = sqrt(data.lambda(1))*ones(size(data.P));
                sdP = sqrt(data.lambda(2))*ones(size(data.P));
            else
                sdY = sqrt(data.lambda(1:data.N))'*ones(1,size(data.P,2));
                sdP = sqrt(data.lambda(data.N+1:end))'*ones(1,size(data.P,2));
            end

            if nargout == 1
                varargout{1} = sdY;
            elseif nargout == 2
                varargout{1} = sdY;
                varargout{2} = sdP;
            end
        end

        function newdata = scaleSNR(data,net,SNR)
        % scales the noise variance to achieve desired SNR. A new dataset is
        % created based on the old data and the new SNR.
        %
        % == Usage ==
        % newdata = scaleSNR(data,net,SNR)
        %
            newdata = datastruct.Dataset(data,net);
            sY = svd(true_response(newdata));
            sE = svd(newdata.E);
            scale = 1/SNR*min(sY)/max(sE);
            newdata.lambda = scale^2*newdata.lambda;
            newdata.E = scale*newdata.E;
        end

        function newdata = variable_std_normalization(data)
        % Standard normalize dataset expression matrix over variables (rows).
        %
        % == Usage ==
        % newdata = variable_std_normalization(data)
        %
            newdata = datastruct.Dataset(data);
            if sum(newdata.E(:)) ~= 0
                SNR = min(svd(true_response(newdata)))/max(svd(newdata.E));
            else
                SNR = 1;
            end
            Y = true_response(newdata);
            m = size(Y, 2);

            mu = mean(Y, 2);
            sigma = std(Y, 1, 2); % population standard deviation
            Mu = repmat(mu, 1, m);
            Sigma = repmat(sigma, 1, m);
            Yhat = (Y - Mu) ./ Sigma;

            sY = svd(Yhat);
            if sum(newdata.E(:)) ~= 0
                sE = svd(newdata.E);
                scale = 1/SNR*min(sY)/max(sE);
            else
                scale = 1;
            end
            newdata.lambda = scale^2*newdata.lambda;
            newdata.E = scale*newdata.E;
            newdata.Y = Yhat + newdata.E;
            newdata.cvP = [];
            newdata.cvY = [];
            [sdY,sdP] = newdata.std();
            newdata.setsdY(sdY);
            newdata.setsdP(sdP);
        end

        function newdata = sample_std_normalization(data)
        % Standard normalize dataset expression matrix over samples (columns).
        %
        % == Usage ==
        % newdata = variable_std_normalization(data)
        %
            newdata = datastruct.Dataset(data);
            if sum(newdata.E(:)) ~= 0
                SNR = min(svd(true_response(newdata)))/max(svd(newdata.E));
            else
                SNR = 1;
            end
            Y = true_response(newdata);
            n = size(Y, 1);

            mu = mean(Y, 1);
            sigma = std(Y, 1, 1); % sample standard deviation
            Mu = repmat(mu, n, 1);
            Sigma = repmat(sigma, n, 1);
            Yhat = (Y - Mu) ./ Sigma;

            sY = svd(Yhat);
            if sum(newdata.E(:)) ~= 0
                sE = svd(newdata.E);
                scale = 1/SNR*min(sY)/max(sE);
            else
                scale = 1;
            end
            newdata.lambda = scale^2*newdata.lambda;
            newdata.E = scale*newdata.E;
            newdata.Y = Yhat + newdata.E;
            newdata.cvP = [];
            newdata.cvY = [];
            [sdY,sdP] = newdata.std();
            newdata.setsdY(sdY);
            newdata.setsdP(sdP);
        end


        function [E,F] = gaussian(data)
        % Generate new gaussian noise matrices E and F with variance lambda for
        % response and/or perturbations.
        %
        % == Usage ==
        % [E,F] = gaussian(data)
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
        % Gives the networks response + noise to input.

            [n,m] = size(data.P);
            if length(varargin) == 1
                if isa(varargin{1},'datastruct.Network')
                    net = varargin{1};
                    X = net.G*(data.P);
                    Y = X + data.E(:,1:m);
                else
                    Y = data.Y;
                end
            else
                Y = data.Y;
            end
        end

        function X = true_response(data,varargin)
        % Gives the networks response to input from data.

            [n,m] = size(data.P);
            if length(varargin) == 1
                if isa(varargin{1},'datastruct.Network')
                    net = varargin{1};
                    X = net.G*(data.P+data.F(:,1:m));
                else
                    X = data.Y - data.E(:,1:m);
                end
            else
                X = data.Y - data.E(:,1:m);
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

        function varargout = eta(data,varargin)
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
            P = data.P;

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

            varargout{1} = etay;
            varargout{2} = etau;

            % data.etay = etay;
            % data.etau = etau;
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
            P = data.P;
            etay = [];
            etau = [];
            for i=1:data.M
                tmpdata = without(data,i);

                [yiU, yiS, yiV] = svd(response(tmpdata,net));
                [uiU, uiS, uiV] = svd(tmpdata.P);
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
            varargout{1} = etay;
            varargout{2} = etau;
        end

        function eta = etay(data)
            [eta,junk] = data.eta();
        end

        function eta = etau(data)
            [junk,eta] = data.eta();
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
                error('Needs to be a struct\n datastruct.Dataset\n datastruct.Experiment or datastruct.Network class')
            end

            inputnames = fieldnames(input);
            names = fieldnames(data);
            for name = inputnames'
                if any(strcmp(name,names))
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

            if nargout == 0
                load@datastruct.Exchange(varargin{:});
                return
            elseif nargout == 1
                varargout{1} = load@datastruct.Exchange(varargin{:});
                return
            elseif nargout == 2
                [varargout{1},varargout{2}] = load@datastruct.Exchange(varargin{:});
                return
            end
        end

        function varargout = fetch(varargin)
        % method to get GeneSPIDER datasets from https://bitbucket.org/sonnhammergrni/gs-datasets
        % currently only supports json format.
        %
        % Several different callings are possible, if URL is a path in the repository (ending in '/')
        % a list of files or directories will be returned.
        %
        % data = datastruct.Dataset.fetch('URL');
        % where URL is the complete path to the file or directory
        % or the extended path beginning with the above URL.
        %
        % data = datastruct.Dataset.fetch('name');
        % where name is the name in the bitbucket repository
        %
        % data = datastruct.Dataset.fetch('option1',value1);
        % set name value pairs to be able to parse options.


            options(1).directurl = '';
            options(1).baseurl = 'https://bitbucket.org/api/1.0/repositories/sonnhammergrni/gs-datasets/raw/';
            options(1).version = 'master';
            options(1).N = 10;
            options(1).name = 'Nordling-ID1446937-D20150825-N10-E15-SNR3291-IDY15968';
            options(1).filetype = '.json';

            if nargin == 0
                if nargout == 0
                    default_url = fullfile(options.baseurl,options.version,['N',num2str(options.N)],'/');
                    fetch@datastruct.Exchange(options,default_url);
                    return
                else
                    default_file = 'Nordling-ID1446937-D20150825-N10-E15-SNR3291-IDY15968.json';
                    obj_data = fetch@datastruct.Exchange(options,default_file);
                end
            else
                obj_data = fetch@datastruct.Exchange(options,varargin{:});
            end

            if isa(obj_data,'cell')
                varargout{1} = obj_data;
                return
            else
                data = datastruct.Dataset();
                data.populate(obj_data);
            end

            if nargout == 1
                varargout{1} = data;
                return
            end

            if nargout == 2
                varargout{1} = data;
                varargout{2} = obj_data;
                return
            end
        end
    end
end
