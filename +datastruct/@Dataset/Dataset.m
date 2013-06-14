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
        SNR       % Signal to noise ratio
        etay      % Sample linear independance from data
        etau      % Perturbation linear independance from data
        info      % Level of informativeness [0,1]
    end

    properties (Hidden = true)
        M           % # variables in A
        N           % # experiments
        created = struct('creator','','time',now,'id','','nexp','')
        tol = eps;
    end

    methods
        function data = Dataset(varargin)
            warning('off','MATLAB:structOnObject')
            if nargin > 0
                for i=1:nargin
                    if isa(varargin{i},'GeneSpider.Network')
                        data.M = size(varargin{i}.A,1);
                        populate(data,varargin{i});
                    elseif isa(varargin{i},'GeneSpider.Experiment')
                        experiment = struct(varargin{i});
                        experiment.Y = trueY(varargin{i});
                        populate(data,experiment);
                        if ispc; data.created.creator = getenv('USERNAME');
                        else; data.created.creator = getenv('USER');
                        end
                        setname(data)
                    elseif isa(varargin{i},'GeneSpider.Dataset')
                        newdata = struct(varargin{i});
                        populate(data,newdata);
                        if ispc; data.created.creator = getenv('USERNAME');
                        else; data.created.creator = getenv('USER');
                        end
                        setname(data)                        
                    elseif isa(varargin{i},'struct')
                        input = varargin{i};
                        populate(data,input);
                    end
                end
                setname(data)
            end
        end

        function SNR = get.SNR(data)
            SNR = min(svd(data.Y))/max(svd(data.E));
        end

        function N = get.N(data)
            N = size(data.P,2);
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
            data.dataset = [namer.creator,'-ID',data.network(regexpi(data.network,'-ID')+3:end),'-D',datestr(namer.time,'yyyymmdd'),'-E',num2str(size(data.P,2)),'-SNR',strrep(num2str(data.SNR),'.','-')];
        end

        function varargout = std(data)
            if numel(data.lambda) == 2,
                sdY = data.lambda(1)*ones(size(data.P));
                sdP = data.lambda(2)*ones(size(data.P));
            else
                sdY = data.lambda(1:data.M)'*ones(1,size(data.P,2));
                sdP = data.lambda(data.M+1:end)'*ones(1,size(data.P,2));
            end
            if nargout == 1
               varargout{1} = sdY
            elseif nargout == 2
               varargout{1} = sdY
               varargout{2} = sdP
            end
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

        function setstdY(data,sdY)
            data.sdY = sdY;
        end

        function setstdP(data,sdP)
            data.sdP = sdP;
        end

        function varargout = cov(data)
        % returns default covariance for the data.
            if numel(data.lambda) == 1
                cvY = data.lambda*eye(data.M);
                cvP = 0*eye(data.M);
            elseif numel(data.lambda) == 2
                cvY = data.lambda(1)*eye(data.M);
                cvP = data.lambda(2)*eye(data.M);
            elseif numel(data.lambda) == data.M
                cvY = diag(data.lambda);
                cvP = 0*eye(size(data.M,1));
            elseif numel(data.lambda) == 2*data.M
                cvY = diag(data.lambda(1:data.M));
                cvP = diag(data.lambda(data.M+1:end));
            end
            if nargout == 1
               varargout{1} = cvY
            elseif nargout == 2
               varargout{1} = cvY
               varargout{2} = cvP
            end
        end

        function setcovY(data,cvY)
           data.cvY = cvY;
        end

        function setcovP(data,cvP)
            data.cvP = cvP;
        end

        function varargout = informativeness(data,net)
        % Gives the level of informativeness of the data set given the network
        % it was produced from.
        %
        % info = informativeness(data,net)
        %
            lambda = data.lambda;
            if numel(lambda) == 2
                o = ones(size(data.P));
                [conf,infotopo] = tools.RInorm(responce(data,net)',data.P',diag(lambda(1:length(lambda)/2))*o',diag(lambda(length(lambda)/2+1:end))*o'+eps);
            else
                o = ones(size(data.P),2);
                [conf,infotopo] = tools.RInorm(responce(data,net)',data.P',(diag(lambda(1:length(lambda)/2))'*o)',(diag(lambda(length(lambda)/2+1:end))'*o)'+eps);
            end

            if nargout == 0
                data.info = sum(sum(logical(net) & infotopo))/sum(sum(logical(net)));
            elseif nargout == 1
                varargout{1} = sum(sum(logical(net) & infotopo))/sum(sum(logical(net)));
            end
        end

        function Y = responce(data,net)
        % Gives the networks responce to input from data.
            if ~isa(net,'GeneSpider.Network')
                error('Must give a GeneSpider.Network to calculate responce')
            end
            Y = net.G*(data.P-data.F) + data.E;
        end

        function newdata = without(data,net,i)
        % creates a Dataset without sample i
        % 
        % == Usage ==
        % newdata = without(data,net,i)
        %           where net is a GeneSpider.Network, and 'i' is the sample
        %           that should be removed.
            N = data.N;
            
            newdata = GeneSpider.Dataset(data,net);

            tmp(1).Y = newdata.Y(:,[1:(i-1) (i+1):N]);
            tmp(1).P = newdata.P(:,[1:(i-1) (i+1):N]);
            tmp(1).E = newdata.E(:,[1:(i-1) (i+1):N]);
            tmp(1).F = newdata.F(:,[1:(i-1) (i+1):N]);
            
            sdY = newdata.sdY;
            if ~isempty(sdY)
                tmp(1).sdY = newdata.sdY;    
            end

            sdP = newdata.sdP;
            if ~isempty(sdP)
                tmp(1).sdP = newdata.sdP;
            end
            
            newdata = populate(newdata,tmp);
        end

        function eta(data,net)
        % calculates eta values for Y and P.
        %
        % eta(data,net)
        %
            Y = responce(data,net);
            P = data.P-data.F;
            etay = [];
            etau = [];
            for i=1:data.N
                tmpdata = without(data,net,i);
                % [Yi,Pi] = without(data,net,i);

                [yiU, yiS, yiV] = svd(tmpdata.Yi);
                [uiU, uiS, uiV] = svd(tmpdata.Pi);
                dyiS = diag(yiS);
                duiS = diag(uiS);

                if length(dyiS) < size(Yi,1)
                    dyiS = [ dyiS; zeros( (size(tmpdata.Yi,1)-length(dyiS)),1 ) ];
                    duiS = [ duiS; zeros( (size(tmpdata.Yi,1)-length(duiS)),1 ) ];
                end

                yv = Y(:, i) / norm( Y(:, i) );
                uv = P(:, i) / norm( P(:, i) );

                yg = dyiS' * abs(yiU' * yv) / sum(dyiS); %yiS(1);
                etay(i) = yg;
                ug = duiS' * abs(uiU' * uv) / sum(duiS); %uiS(1);
                etau(i) = ug;
            end
            data.etay;
            data.etau;
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