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
            net.interrampatteness = cond(full(A));
            net.NetworkComponents = graphconncomp(sparse(sign(A)), 'Directed', true);
            net.tauG = min(sort(1./abs(real(eig(-pinv(full(A)))))));

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
            s = svd(full(net.A));
        end

        function Y = response(net,data)
        % Gives the network response to input from data.
            if ~isa(data,'datastruct.Dataset')
                error('Must give a datastruct.Dataset to calculate response')
            end
            Y = net.G*(data.P-data.F) + data.E;
        end

        function net = populate(net,input)
        % populate the Network object with fields of the network struct.
        %
            if ~isa(input,'struct') && ~isa(input,'datastruct.Network')
                error('Needs to be a struct or datastruct.Network object')
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
    end
end