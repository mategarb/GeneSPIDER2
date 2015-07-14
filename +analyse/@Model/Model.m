classdef Model

    properties (SetAccess = private)
        network = '';     % identifier for network used
        interrampatteness % interrampatteness
        NetworkComponents % Number of strong components in network
        AvgPathLength     % Average path length
        tauG              % Time Constant
        
        % N                 % Number of variables in model

        %% top_entropy       % Structural/topological entropy
    end

    properties (SetAccess = public)
        tol = eps;         % tolerance level
        type = 'directed'; % type of network
        alpah = 0.01;      % confidence level
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
            % analysis.N                 = analyse.Model.size(mat);
        end
    end

    methods (Static)
        function A = give_matrix(model,varargin)
            if isa(model,'datastruct.Network') && nargin < 2
                A = model.A;
            elseif isa(model,'datastruct.Network') && strcmpi(varargin{1},'inv')
                A = model.G;
            else
                A = model;
            end
        end

        function network = identifier(model)
            if isa(model,'datastruct.Network')
                network = model.network;
            else
                network = '';
            end
        end

        % function N = size(model)
        %     A = analyse.Model.give_matrix(model);
        %     N = size(A,1);
        % end

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

        function NetworkComponents = graphconncomp(model)
            A = analyse.Model.give_matrix(model);
            NetworkComponents = graphconncomp(sparse(sign(A)), 'Directed', true);
        end

        function tauG = time_constant(model)
            G = analyse.Model.give_matrix(model,'inv');
            tauG = min(sort(1./abs(real(eig(G)))));
        end
    end
end
