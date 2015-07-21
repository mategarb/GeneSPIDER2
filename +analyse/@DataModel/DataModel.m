classdef DataModel

    methods
        function analysis = DataModel()
            
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
    end
    
end