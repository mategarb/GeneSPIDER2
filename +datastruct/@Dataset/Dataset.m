classdef Dataset < hgsetget
% The Dataset class is used to store a GeneSpider.Dataset
% complementary to a GeneSpider.Network
%
%   Input arguments: Dataset([Network ,Experiment])
%   ================
%   (none) :      Initiate Dataset
%   Network :     A Network object.
%   Experiment :  An Experiment object or a structure with named fields.
%                 If an Experiment is not supplied
%                 then data needs to be added manualy later with the function populate
%
%   ================ functionality
%

    properties (SetAccess = private)
        dataset   % Name of the data set
        network   % Name of complementary network
        P         % True perturbations
        F         % Perturbation noise
        cvP = []; % Covariance of P
        sdP = []; % measurement point variation of P
        Y         % True expression response
        E         % Expression noise
        cvY = []; % Covariance of noisy Y
        sdY = []; % measurement point variation of Y
        lambda    % Noise variance
        SNR       % Signal to noise ratio, \sigma_N(Y)/\sigma_1(E)
        SNRnu     % Signal to noise ratio, \argmin_i min(norm(\Fi_i)/norm(\nu_i))
        SNRm      % Signal to noise ratio,
        SNRv      % Signal to noise ratio,
        etay      % Linear independance from data per sample.
        etau      % Linear independance from data per perturbation.
        info      % Level of informativeness [0,1]
    end

    properties (Hidden = true)
        N           % # variables in A
        M           % # experiments
        created = struct('creator','','time',now,'id','','nexp','');
        tol = eps;
        alpha = 0.05; % Confidence
    end

    methods
        function data = Dataset(varargin)
            warning('off','MATLAB:structOnObject')
            if nargin > 0
                for i=1:nargin
                    if isa(varargin{i},'GeneSpider.Network')
                        populate(data,varargin{i});
                    elseif isa(varargin{i},'GeneSpider.Experiment')
                        experiment = struct(varargin{i});
                        experiment.Y = trueY(varargin{i});
                        populate(data,experiment);
                        if ispc; data.created.creator = getenv('USERNAME');
                        else; data.created.creator = getenv('USER');
                        end
                        data.created.id = num2str(round(cond(data.Y)*10000));
                    elseif isa(varargin{i},'GeneSpider.Dataset')
                        newdata = struct(varargin{i});
                        populate(data,newdata);
                        if ispc; data.created.creator = getenv('USERNAME');
                        else; data.created.creator = getenv('USER');
                        end
                        data.created.id = num2str(round(cond(data.Y)*10000));
                    elseif isa(varargin{i},'struct')
                        input = varargin{i};
                        populate(data,input);
                    end
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

        function SNR = get.SNR(data)
            SNR = min(svd(data.Y))/max(svd(data.E));
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

        function SNR = get.SNRnu(data)
            snr = [];
            for i=1:data.N
                snr(i) = norm(data.Y(i,:))/norm(data.E(i,:));
            end
            SNR = min(snr);
        end

        function SNRm = get.SNRm(data)
            alpha = data.alpha;
            sigma = min(svd(data.Y'));
            SNRm = sigma/sqrt(chi2inv(1-alpha,prod(size(data.P)))*data.lambda(1));
        end

        function SNRv = get.SNRv(data)
            alpha = data.alpha;
            for i=1:data.N
                snr(i) = norm(data.Y(i,:))/sqrt(chi2inv(1-alpha,prod(size(data.P)))*data.lambda(1));
            end
            SNRv = min(snr);
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
            data.dataset = [namer.creator,'-ID',data.network(regexpi(data.network,'-ID')+3:end),'-D',datestr(namer.time,'yyyymmdd'),'-E',num2str(size(data.P,2)),'-SNR',num2str(round(data.SNRm*1000)),'-IDY',namer.id];
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
            s = svd(data.Y);
            % preLambda = data.lambda(1);
            lambda = min(s)^2/(chi2inv(1-data.alpha,prod(size(data.P)))*SNRm^2);
            data.lambda = lambda;
            % data.E = sqrt(lambda/preLambda).*data.E;
        end

        function newdata = scaleSNR(data,net,SNR)
        % scales the noise variance to achieve desired SNR. A new dataset is
        % created based on the old data and the new SNR.
        %
        % == Usage ==
        % newdata = scaleSNR(data,net,SNR)
        %
            newdata = GeneSpider.Dataset(data,net);
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

        function varargout = informativeness(data,varargin)
        % Gives the level of informativeness of the data set given the network
        % it was produced from.
        %
        % [info{, uninfo}] = informativeness(data,net)
        %
            net = [];
            if length(varargin) == 1
                if isa(varargin{1},'GeneSpider.Network')
                    net = varargin{1};
                end
            end
            lambda = data.lambda;
            if numel(lambda) == 2
                o = ones(size(data.P));
                [conf,infotopo] = tools.RInorm(response(data,net)',data.P',diag(lambda(1:length(lambda)/2))*o',diag(lambda(length(lambda)/2+1:end))*o'+eps,data.alpha);
            else
                o = ones(size(data.P),2);
                [conf,infotopo] = tools.RInorm(response(data,net)',data.P',(diag(lambda(1:length(lambda)/2))'*o)',(diag(lambda(length(lambda)/2+1:end))'*o)'+eps,data.alpha);
            end

            if nargout == 0
                data.info = sum(sum(logical(net) & infotopo))/sum(sum(logical(net)));
                return
            elseif nargout >= 1
                varargout{1} = sum(sum(logical(net) & infotopo))/sum(sum(logical(net)));
            end
            if nargout >= 2
                varargout{2} = sum(conf(~logical(net)))/sum(sum(~logical(net)));
            end
            if nargout >= 3
                varargout{3} = conf;
            end
        end

        function Y = response(data,varargin)
        % Gives the networks response to input from data.

            [n,m] = size(data.P);
            if length(varargin) == 1
                if isa(varargin{1},'GeneSpider.Network')
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
        %           and net is a GeneSpider.Network.
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
                newdata = GeneSpider.Dataset(data,net);
            else
                newdata = GeneSpider.Dataset(data);
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
                if isa(varargin{1},'GeneSpider.Network')
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
                if isa(varargin{1},'GeneSpider.Network')
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
        % based on eta limit, returns samples ok to inkclude in LOCO
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
                included = [and(data.etay >= SY(end), data.etau >= SP(end))];
            elseif exist('etaLim','var')
                if numel(etaLim) == 1
                    included = [and(data.etay >= etaLim, data.etau >= etaLim)];
                else
                    included = [and(data.etay >= etaLim(1), data.etau >= etaLim(2))];
                end
            end
        end

        function varargout = populate(data,input)
        % populate the Dataset object with matching fields of the input.
        %
        % == Usage ==
        % {data =} populate(data,input)
        %          where input can be a struct, GeneSpider.Dataset,
        %          GeneSpider.Experiment or GeneSpider.Network

            if ~isa(input,'struct') && ~isa(input,'GeneSpider.Dataset') && ~isa(input,'GeneSpider.Experiment') && ~isa(input,'GeneSpider.Network')
                error('Needs to be a struct,\n GeneSpider.Dataset,\n GeneSpider.Experiment or GeneSpider.Network class')
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
        % saves a GeneSpider.Dataset to file, either mat or xml.
        % The name of the file will be the name of the data set.
        %
        %   Input Arguments: savedata(dataset[,path,<fileext>])
        %   ================
        %   (none) :        Save data file to current directory as a .mat file.
        %   path :          Path to directory.
        %   fileext :       File extension: .xml or (.mat)
        %
            warning('off','MATLAB:structOnObject')
            fending = '.mat';
            savepath = './';
            dataset = struct(data);
            if nargin > 1
                if ~isa(varargin{1},'char')
                    error('Input arguments must be a string')
                end
                savepath = varargin{1};
            end

            if nargin == 3
                if ~isa(varargin{2},'char')
                    error('Input arguments must be a string')
                end
                fending = varargin{2};
                if strcmp(fending,'.xml')
                    if exist('mat2xml') ~= 2
                        error('Save method for xml files does not seem to exist')
                    end
                end
            end

            name = dataset.dataset;
            savevar = 'dataset';

            if strcmp(fending,'.xml')
                xmlString = simplify_mbml( spcharout( mat2xml(dataset,savevar)) );
                xmlwrite(fullfile(savepath,[name,'.xml']), str2DOMnode(xmlString));
            elseif strcmp(fending,'.mat')
                save(fullfile(savepath,name),savevar)
            else
                error('unknown file extension')
            end
        end
    end

    methods (Static)
        function varargout = load(varargin)
        % Load a dataset file back in to a Dataset object
        % dataset = GeneSpider.Dataset.loaddata(['path/file'] or [path,file]);
        %
        %   Input Arguments: GeneSpider.Dataset.loaddata([path,file])
        %   ================
        %   (none) :        Outputs a list of datasets availible in the current directory.
        %   path :          Path to directory or full path with filename. If no file is specified
        %                   it will output a list of availible data sets in that directory.
        %   file :          Filename or number of its place in the list of files.
        %                   Have to include the path input variable.
        %
        %   Output Arguments: net / list
        %   ================
        %   net :           Populate the Network object with the loaded file.
        %   list :          If no file is specified a list of availible datasets
        %                   is returned for the given directory.
        %

            lpath = pwd;
            lfile = [];
            if nargin == 1
                if isa(varargin{1},'double')
                    lfile = varargin{1};
                else
                    if exist(varargin{1}) == 2
                        [p,f,e] = fileparts(varargin{1});
                        lpath = p;
                        lfile = [f,e];
                    elseif exist(varargin{1}) == 7
                        lpath = varargin{1};
                    else
                        error('Unknown path or file')
                    end
                end
            elseif nargin == 2
                lpath = varargin{1};
                if exist(lpath) ~= 7
                    error('Unknown path')
                end
                if isa(varargin{2},'double')
                    lfile = varargin{2};
                else
                    if exist(fullfile(lpath,varargin{2})) == 2
                        lfile = varargin{2};
                    else
                        error('Unknown file')
                    end
                end
            elseif nargin > 2
                error('wrong number of input arguments.')
            end

            if isa(lfile,'double')
                datasets = dir(lpath);
                j = 0;
                for i = 1:length(datasets)
                    if ~datasets(i).isdir
                        [pa,fi,fext] = fileparts(datasets(i).name);
                        if strcmp(fext,'.mat') || strcmp(fext,'.xml')
                            j = j+1;
                            output{j} = datasets(i).name;
                        end
                    end
                end
                if nargout == 1 && isempty(lfile)
                    varargout{1} = output;
                    return
                elseif isempty(lfile)
                    for j=1:length(output)
                        fprintf('%d %s\n',j,output{j});
                    end
                    return
                else
                    lfile = output{lfile};
                end
            end

            data = GeneSpider.Dataset;
            fetchfile = fullfile(lpath,lfile);
            [p,f,e] = fileparts(fetchfile);
            if strcmp(e,'.mat')
                load(fetchfile);
            elseif strcmp(e,'.xml')
                [MAT,dataset] = xml2mat(fetchfile);
                eval([dataset,'=MAT;']);
            end

            populate(data,dataset);
            if nargout == 1
                varargout{1} = data;
            elseif nargout == 2
                varargout{1} = data;
                varargout{2} = dataset;
            end
        end
    end
end