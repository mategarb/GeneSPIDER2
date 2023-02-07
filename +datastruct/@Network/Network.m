classdef Network < datastruct.Exchange
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
        network           % Uniq name of network
        A                 % network
        G                 % Static gain model
    end

    properties
        names = {};
        description = '';
    end

    properties (Hidden = true)
        created = struct('creator','', 'time',now,'id','','nodes','','type','unknown','sparsity','')
        tol = eps;
        N                 % # Nodes in network
    end

    methods
        function net = Network(varargin)

            net = net@datastruct.Exchange();
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
            net.G = -pinv(full(A));

            net.created.id = num2str(round(cond(full(A))*10000));
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

        function show(net)
            networkProperties = {
                'Name'                 net.network;
                'Description'          net.description;
                'Sparseness'           nnz(net)/numel(net)
                '# Nodes'              size(net.A,1)
                '# links'              nnz(net)};

            f = figure();
            t1 = uitable('Parent', f, 'Position', [10 10 730 310],'ColumnWidth',{65});
            set(t1, 'Data', net.A,'ColumnName',net.names,'RowName',net.names,'BackgroundColor',[1,1,1]);
            f = figure();
            t2 = uitable('Parent', f,'ColumnWidth',{200,480});
            set(t2, 'Data', networkProperties,'BackgroundColor',[1,1,1],'ColumnName',{'Property','Value'});
        end

        function varargout = view(net)
        % make a rough graphical network plot with biograph
            Ap = net.A;
            Ap(logical(eye(net.N))) = 0;
            h = view(biograph(Ap,net.names,'ShowWeights','on'));

            if nargout == 1
                varargout{1} = h;
            end
        end

        function SA = sign(net)
            SA = sign(net.A);
        end

        function LA = logical(net)
            LA = logical(net.A);
        end

        function sA = size(net,varargin)
            if ~isempty(varargin)
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
            s = svd(full(net.A));
        end

        function Y = mtimes(net,p)
        % generate a mapping of a perturbation through the network returning SS response
        %
        % y = net*p;
        %

            if isrow(p)
                p = p';
            end
            Y = net.G*p;

        end

        function net = populate(net,input)
        % populate the Network object with fields of the network struct.
        %
            if ~isa(input,'struct') && ~isa(input,'datastruct.Network')
                error('Needs to be a struct or datastruct.Network object')
            end
            inputnames = fieldnames(input);
            allnames = fieldnames(net);
            for name = inputnames'
                if any(find(strcmp(name,allnames)))
                    net.(name{1}) = input.(name{1});
                end
            end
        end

        function save(net,varargin)
            save@datastruct.Exchange(net,varargin{:});
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
        % method to get GeneSPIDER networks from https://bitbucket.org/sonnhammergrni/gs-networks
        %
        % Several different callings are possible, if URL is a path in the repository (ending in '/')
        % a list of files or directories will be returned.
        %
        % net = datastruct.Network.fetch('URL');
        % where URL is the complete path to the file or directory
        % or the extended path beginning with the above URL.
        %
        % net = datastruct.Network.fetch('name');
        % where name is the name in the bitbucket repository.
        %
        % net = datastruct.Network.fetch('option1',value1);
        % set name value pairs to be able to parse options.

            options(1).directurl = '';
            options(1).baseurl = 'https://bitbucket.org/sonnhammergrni/gs-networks/raw/';
            options(1).version = 'master';
            options(1).type = 'random';
            options(1).N = 10;
            options(1).name = '';
            options(1).filelist = false;
            options(1).filetype = '';

            if nargin == 0
                if nargout == 0
                    default_url = fullfile(options.baseurl,options.version,options.type,['N',num2str(options.N)],'/');
                    fetch@datastruct.Exchange(options,default_url);
                    return
                else
                    default_file = 'Nordling-D20100302-random-N10-L25-ID1446937.json';
                    obj_data = fetch@datastruct.Exchange(options,default_file);
                end
            else
                if nargout == 0
                    fetch@datastruct.Exchange(options,varargin{:});
                    return
                else
                    obj_data = fetch@datastruct.Exchange(options,varargin{:});
                end
            end

            if isa(obj_data,'cell')
                varargout{1} = obj_data;
                return
            else
                net = datastruct.Network();
                net.populate(obj_data);
            end

            if nargout == 1
                varargout{1} = net;
                return
            elseif nargout == 2
                varargout{1} = net;
                varargout{2} = obj_data;
                return
            end
        end
    end
end
