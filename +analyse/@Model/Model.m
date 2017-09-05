classdef Model < analyse.DataModel
% Model calculates properties related to the supplied network.
%
%
%
    properties (SetAccess = private)
        network = '';     % identifier for network used
        interampatteness  % interampatteness
        proximity_ratio   % the small world tendency
        networkComponents % Number of strong components in network
        medianPathLength  % Median path length
        meanPathLength    % Mean path length
        tauG              % Time Constant
        CC                % average clustering coefficient
        DD                % average degree distribution

        %% top_entropy       % Structural/topological entropy
    end

    methods
        function analysis = Model(mat,varargin)

            analysis = analysis@analyse.DataModel();
            if nargin > 1
                ind = find(strcmpi('tol',varargin));
                if ind
                    analysis.tol = varargin{ind+1};
                end
            end

            analysis = analyse_model(analysis,mat,varargin);

        end

        function analysis = analyse_model(analysis,mat,varargin);
            analysis.network                                    = analyse.Model.identifier(mat);
            analysis.interampatteness                           = analyse.Model.cond(mat);
            analysis.networkComponents                          = analyse.Model.graphconncomp(mat);
            [analysis.medianPathLength,analysis.meanPathLength] = analyse.Model.median_path_length(mat);
            analysis.tauG                                       = analyse.Model.time_constant(mat);
            analysis.CC                                         = nanmean(analyse.Model.clustering_coefficient(mat));
            analysis.DD                                         = mean(analyse.Model.degree_distribution(mat));
            analysis.proximity_ratio                            = analyse.Model.calc_proximity_ratio(mat);
        end
    end

    methods (Static)
        function varargout = alpha(varargin)
        % significance level (default = 0.01)
            val = alpha@analyse.DataModel(varargin{:});
            if nargout == 1
                varargout{1} = val;
                return
            end
            if nargin == 0
                disp(val)
            end
        end

        function varargout = type(varargin)
        % 'directed' or 'undirected'
            val = type@analyse.DataModel(varargin{:});
            if nargout == 1
                varargout{1} = val;
                return
            end
            if nargin == 0
                disp(val)
            end
        end

        function varargout = tol(varargin)
        % if a tolernace value for computations is needed it can be set here.
            val = tol@analyse.DataModel(varargin{:});
            if nargout == 1
                varargout{1} = val;
                return
            end
            if nargin == 0
                disp(val)
            end
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

        function interampatteness = cond(model)
            A = analyse.Model.give_matrix(model);
            interampatteness = cond(full(A));
        end

        function varargout = median_path_length(model,varargin)
        % Calculate mean and median path length
        % [median_pl,mean_pl,median_pl,mean_pl] = analyse.Model.median_path_length(network)
        %
        % Any additional arguments are passed directly to matlabs
        % graphallshortestpaths function
        %

            A = analyse.Model.give_matrix(model);
            pl = graphallshortestpaths(logical(sparse(A)),'directed',strcmpi(analyse.Model.type,'directed'),varargin{:});
            pl = pl(~eye(size(A)));
            pl = pl(~isinf(pl));

            if nargout > 0
                mpl = median(pl);
                varargout{1} = mpl;
            end
            if nargout > 1
                apl = mean(pl);
                varargout{2} = apl;
            end
            
            if nargout > 2
                pl = graphallshortestpaths(sparse(A),'directed',strcmpi(analyse.Model.type,'directed'),varargin{:});
                pl = pl(~eye(size(A)));
                pl = pl(~isinf(pl));
                
                mpl = median(pl);
                varargout{3} = mpl;
            end
            if nargout > 3
                apl = mean(pl);
                varargout{4} = apl;
            end
            
            if nargout > 4
                numinf = length(isinf(pl));
                varargout{5} = numinf;
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
        % Calculates minimum time constant of the model G
        % tauG = time_constant(G)
        %
        % G is a datastruct.Network object.
        %
            G = analyse.Model.give_matrix(model,'inv');
            tauG = min(sort(1./abs(real(eig(G)))));
        end

        function varargout = clustering_coefficient(model)
        % Calculate clustering coefficient
        %
        % [c_out,c_in] = analyse.Model.clustering_coefficient(net)
        %
        % c_out is the clustering coefficient considering outgoing links of the node in question.
        % c_in  is the clustering coefficient considering incoming links of the node in question.
        %
        % For undirected networks only one clustering coeficient is calculated.
        %
        % net is the network or matrix to be tested.
        %
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

                end
                C_out(find(isnan(C_out))) = 0;
                C_in(find(isnan(C_in))) = 0;
                varargout{1} = C_out;
                varargout{2} = C_in;
            elseif strcmpi(analyse.Model.type,'undirected')
                A = logical(A + A');
                for i=1:size(A,1)
                    local_vertices = find(A(:,i));
                    subnet = A(local_vertices,local_vertices);
                    C(i) = nnz(subnet)/(prod(size(subnet))-size(subnet,1));
                end
                varargout{1} = C;
            end
        end

        function varargout = degree_distribution(model)
        % Calculate degree distribution
        %
        % [out_deg,in_deg] = analyse.Model.degree_distribution(net)
        %
        % out_deg is the degree distribution of outgoing links.
        % in_deg  is the degree distribution of incoming links.
        %
        % For undirected networks only one clustering coeficient is calculated.
        %
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
                varargout{2} = degree_dist;
            end
        end

        function varargout = calc_proximity_ratio(model)
        % Calculates the proximity ratio or "smallworldness" tendency of the network
        % > 1 means higher than for a random network
        % < 1 means lower than for a random network
        %
            A = analyse.Model.give_matrix(model);
            A(find(eye(size(A)))) = 0; % remove self loops
            A = logical(A);
            N = size(A,1);

            [Lt,L] = analyse.Model.median_path_length(A);
            md = nnz(A)/N;
            Lr = log(N)/log(md);

            C = mean(analyse.Model.clustering_coefficient(A));
            Cr = nnz(A)/(N*(N-1));

            smallworldness = (C/Cr)*(Lr/L);

            varargout{1} = smallworldness;
        end
    end
end
