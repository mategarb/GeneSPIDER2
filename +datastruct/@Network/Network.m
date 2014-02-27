classdef Network < hgsetget
% Stores a network matrix, A and calculates the most important network properties.
%
%   net = Network(A);
%
%   Some basic functionality can be done directly on this object:
%   ==============
%   svd(net)
%   sign(net)
%   logical(net)
%   nnz(net)
%   size(net)
%
%   ============== show
%   show() the Network will launch a graphics window and display a matrix
%   of the network and a table with network properties.
%
    properties (SetAccess = private)
        A                 % network
        G                 % Static gain model
        network           % Uniq name of network
        interrampatteness % interrampatteness
        NetworkComponents % Number of strong components in network
        tauG              % Time
    end

    properties
        names = {};
        desc = '';
    end

    properties (Hidden = true)
        created = struct('creator','', 'time',now,'id','','nodes','','type','','sparsity','')
        tol = eps;
        N                 % # Nodes in network
    end

    methods
        function net = Network(varargin)
            if nargin == 2
                setA(net,varargin{1});
                net.created.type = varargin{2};
                setname(net);
            elseif nargin == 1
                setA(net,varargin{1});
                setname(net);
            end
        end

        function setA(net,A)
            net.A = A;
            net.G = -pinv(A);
            net.interrampatteness = cond(A);
            net.NetworkComponents = graphconncomp(sparse(sign(A)), 'Directed', true);
            net.tauG = min(sort(1./abs(real(eig(-pinv(A))))));

            net.created.id = num2str(round(cond(A)*10000));
            net.created.nodes = num2str(size(A,1));
            net.created.sparsity = nnz(net);
            if ispc; net.created.creator = getenv('USERNAME');
            else; net.created.creator = getenv('USER');
            end
        end

        function setname(net,varargin)
            if nargin == 1
                namestruct = net.created;
            elseif nargin == 2
                namestruct = varargin{1};
            end
            if ~isa(namestruct,'struct')
                error('input argument must be name/value pairs in struct form')
            end

            namer = net.created;
            optNames = fieldnames(namer);
            inpNames = fieldnames(namestruct);

            for i=1:length(inpNames)
                if any(strcmp(inpNames{i},optNames))
                    namer.(inpNames{i}) = namestruct.(inpNames{i});
                end
            end

            net.network = [namer.creator,'-D',datestr(namer.time,'yyyymmdd'),'-',namer.type,'-N',namer.nodes,'-L',num2str(nnz(net)),'-ID',namer.id];
        end

        function N = get.N(net)
            N = size(net,1);
        end

        function names = get.names(net)
            names = net.names;
            if isempty(names)
                for i=1:net.N
                    names{i} = sprintf(['G%0',num2str(floor(log10(net.N))+1),'d'],i);
                end
            else
                names = net.names;
            end
        end

        function desc = get.desc(net)
            desc = sprintf([net.desc,'\n']);
        end

        function show(net)
            networkProperties = {
                'Name'                 net.network;
                'Description'          net.desc;
                'interrampatteness'    net.interrampatteness;
                '# Network Components' net.NetworkComponents;
                '\tau G'               net.tauG
                'Sparseness'           nnz(net)/prod(size(net))
                '# Nodes'              size(net.A,1)
                '# links'              nnz(net)};

            f = figure();
            t1 = uitable('Parent', f, 'Position', [10 10 730 310],'ColumnWidth',{65});
            set(t1, 'Data', net.A,'ColumnName',net.names,'RowName',net.names,'BackgroundColor',[1,1,1]);
            f = figure();
            t2 = uitable('Parent', f,'ColumnWidth',{200,480});
            set(t2, 'Data', networkProperties,'BackgroundColor',[1,1,1],'ColumnName',{'Property','Value'});
        end

        function SA = sign(net)
            SA = sign(net.A);
        end

        function LA = logical(net)
            LA = logical(net.A);
        end

        function sA = size(net,varargin)
            if length(varargin) > 0
                sA = size(net.A,varargin{1});
            else
                sA = size(net.A);
            end
        end

        function nnzA = nnz(net)
            nnzA = nnz(net.A);
        end

        function s = svd(net)
        % returns the singular values (only) of A. To get the singular vectors use svd(net.A)
            s = svd(net.A);
        end

        function info = informativeness(net,data)
            o = ones(size(net));
            [conf,infotopo] = tools.RInorm(response(net,data)',data.P',diag(lambda(1:length(lambda)/2))*o,diag(lambda(length(lambda)/2+1:end))*o+eps);

            info = sum(sum(logical(net) & infotopo));
        end

        function Y = response(net,data)
        % Gives the network response to input from data.
            if ~isa(data,'GeneSpider.Dataset')
                error('Must give a GeneSpider.Dataset to calculate response')
            end
            Y = net.G*(data.P-data.F) + data.E;
        end

        function net = populate(net,input)
        % populate the Network object with fields of the network struct.
        %
            if ~isa(input,'struct') && ~isa(input,'GeneSpider.Network')
                error('Needs to be a struct or Genespider.Network class')
            end
            inputnames = fieldnames(input);
            names = fieldnames(net);
            for name = inputnames'
                if any(strmatch(name,names,'exact'))
                    net.(name{1}) = input.(name{1});
                end
            end
        end

        function save(net,varargin)
        % saves a GeneSpider.Network to file, either mat or xml.
        % The name of the file will be the name of the data set.
        %
        %   Input Arguments: savenet(network[,path,<fileext>])
        %   ================
        %   (none):        Save network file to current directory as a .mat file.
        %   path:          Path to directory.
        %   fileext:       File extension: .xml or (.mat)
        %

            warning('off','MATLAB:structOnObject')
            fending = '.mat';
            savepath = './';
            network = struct(net);
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

            name = network.network;
            savevar = 'network';

            if strcmp(fending,'.xml')
                xmlString = simplify_mbml( spcharout( mat2xml(network,savevar)) );
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
        % Load a network data-file back in to a Network object
        % net = GeneSpider.Network.loadnet(['path/file'] or [path,file]);
        %
        %   Input Arguments: GeneSpider.Network.loadnet([,path,file])
        %   ================
        %   (none) :        Outputs a list of datasets availible in the current directory.
        %   path :          Path to direcotry or full path with filename. If no file is specified
        %                   it will output a list of availible data sets in that directory.
        %   file :          Filename or number of its place in the list of files.
        %                   Have to include the path input variable.
        %
        %   Output Arguments:
        %   ================
        %   net :           Populate the Network object with the loaded file.
        %   list :          If no file is specified a list of availible datasets is returned.
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
                networks = dir(lpath);
                j = 0;
                for i = 1:length(networks)
                    if ~networks(i).isdir
                        [pa,fi,fext] = fileparts(networks(i).name);
                        if strcmp(fext,'.mat') || strcmp(fext,'.xml')
                            j = j+1;
                            output{j} = networks(i).name;
                        end
                    end
                end
                if nargout == 1 && isempty(lfile)
                    varargout{1} = output;
                    return
                elseif isempty(lfile)
                    if exist('output','var')
                        for j=1:length(output)
                            fprintf('%d %s\n',j,output{j});
                        end
                        return
                    else
                        error('No mat- or xml-files to list in dir %s',lpath)
                    end
                else
                    lfile = output{lfile};
                end
            end

            net = GeneSpider.Network;
            fetchfile = fullfile(lpath,lfile);
            [p,f,e] = fileparts(fetchfile);
            if strcmp(e,'.mat')
                load(fetchfile);
            elseif strcmp(e,'.xml')
                [MAT,network] = xml2mat(fetchfile);
                eval([network,'=MAT;']);
            end

            populate(net,network);
            if nargout == 1
                varargout{1} = net;
            elseif nargout == 2
                varargout{1} = net;
                varargout{2} = network;
            end
        end
    end
end