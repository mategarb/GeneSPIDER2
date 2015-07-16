classdef Model

    properties (SetAccess = private)
        network = '';     % identifier for network used
        interrampatteness % interrampatteness
        NetworkComponents % Number of strong components in network
        AvgPathLength     % Average path length
        tauG              % Time Constant
        CC                % average clustering coefficient
        DD                % average degree distribution

        %% top_entropy       % Structural/topological entropy
    end

    methods
        function analysis = Model(mat,varargin)

            if nargin > 1
                ind = find(strcmpi('tol',varargin));
                if ind
                    analysis.tol = varargin{ind+1};
                end
            end

            analysis = analyse_model(analysis,mat,varargin);

        end

        function analysis = analyse_model(analysis,mat,varargin);
            analysis.network           = analyse.Model.identifier(mat);
            analysis.interrampatteness = analyse.Model.cond(mat);
            analysis.NetworkComponents = analyse.Model.graphconncomp(mat);
            analysis.AvgPathLength     = analyse.Model.average_path_length(mat);
            analysis.tauG              = analyse.Model.time_constant(mat);
            analysis.CC                = nanmean(analyse.Model.clustering_coefficient(mat));
            analysis.DD                = mean(analyse.Model.degree_distribution(mat));
        end
    end

    methods (Static)
        function val = alpha(newval)
        % significance level (default = 0.01)
            persistent currentval;

            if isempty(currentval)
                currentval = 0.01;
            end
            if nargin >= 1
                currentval = newval;
            end
            val = currentval;
        end

        function val = type(newval)
        % 'directed' or 'undirected'
            persistent currentval;
            if isempty(currentval)
                currentval = 'directed';
            end
            if nargin >= 1
                currentval = newval;
            end
            val = currentval;
        end

        function val = tol(newval)
        % if a tolernace value for computations is needed it can be set here.
            persistent currentval;
            if isempty(currentval)
                currentval = eps;
            end
            if nargin >= 1
                currentval = newval;
            end
            val = currentval;
        end

        function A = give_matrix(model,varargin)
        % make sure to return a matrix of the network
            if isa(model,'datastruct.Network') && nargin < 2
                A = model.A;
            elseif isa(model,'datastruct.Network') && strcmpi(varargin{1},'inv')
                A = model.G;
            else
                A = model;
            end
        end

        function network = identifier(model)
        % name of network used.
            if isa(model,'datastruct.Network')
                network = model.network;
            else
                network = '';
            end
        end

        function interrampatteness = cond(model)
            A = analyse.Model.give_matrix(model);
            interrampatteness = cond(full(A));
        end

        function varargout = average_path_length(model,varargin)
            A = analyse.Model.give_matrix(model);
            pl = graphallshortestpaths(logical(sparse(A)),varargin{:});
            apl = mean(pl(~eye(size(A))));

            if nargout > 0
                varargout{1} = apl;
            end

        end

        function NetworkComponents = graphconncomp(model,varargin)
            A = analyse.Model.give_matrix(model);

            if strcmpi(analyse.Model.type,'directed')
                flag = true;
            else
                flag = false;
            end

            NetworkComponents = graphconncomp(sparse(sign(A)), 'directed', flag);
        end

        function tauG = time_constant(model)
            G = analyse.Model.give_matrix(model,'inv');
            tauG = min(sort(1./abs(real(eig(G)))));
        end

        function varargout = clustering_coefficient(model)
        % Calculate clustering coefficient
            A = analyse.Model.give_matrix(model);
            A(find(eye(size(A)))) = 0; % remove self loops

            if strcmpi(analyse.Model.type,'directed')
                for i=1:size(A,1)
                    local_vertices = find(A(:,i));
                    subnet = A(local_vertices,local_vertices);
                    C_out(i) = nnz(subnet)/(prod(size(subnet))-size(subnet,1));

                    local_vertices = find(A(i,:));
                    subnet = A(local_vertices,local_vertices);
                    C_in(i) = nnz(subnet)/(prod(size(subnet))-size(subnet,1));

                    varargout{1} = C_out;
                    varargout{2} = C_in;
                end
            elseif strcmpi(analyse.Model.type,'undirected')
                A = logical(A + A');
                for i=1:size(A,1)
                    local_vertices = find(A(:,i));
                    subnet = A(local_vertices,local_vertices);
                    C(i) = nnz(subnet)/(prod(size(subnet))-size(subnet,1));
                    varargout{1} = C;
                end
            end
        end

        function varargout = degree_distribution(model)
        % Calculate degree distribution
            A = analyse.Model.give_matrix(model);
            A(find(eye(size(A)))) = 0; % remove self loops
            A = logical(A);

            if strcmpi(analyse.Model.type,'directed')
                out_degree = sum(A,1);
                in_degree = sum(A',1);
                varargout{1} = out_degree;
                varargout{2} = in_degree;
            elseif strcmpi(analyse.Model.type,'undirected')
                A = logical(A + A');
                degree_dist = sum(A);
                varargout{1} = degree_dist;
            end
        end
    end
end
