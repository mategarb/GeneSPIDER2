classdef CompareModels
% CompareModels calculates difference measures between weighted network adjacency matrices.
%
%   Input Arguments: CompareModels(A, Alist, selected)
%   ================
%   (none)     Initiate object for later use.
%   A        = The 'true' weighted adjacency matrix of the network, to which all
%              other network matrices are compared.
%   Alist    = A network matrix or 3d array of network matrices that are
%              compared to the 'true' one.
%
%   Output Arguments:
%   ================
%   obj = CompareModels object.
%
%   Use Examples:
%   ================
%
%   % Create sparse random network:
%   A = full(sprandn(10,10,0.25));
%
%   % Create random Alist networks:
%   Alist = randn(10,10,20);
%   Alist(Alist < 0.75) = 0;
%
%
%   M = analyse.CompareModels();              % initiate comparision
%   % set(M,'tol',value)                      % set tolerance {eps}
%   M.A = A;                                  % set true A
%   <results> = compare(M,Alist[,selected]);  % Update M with comparisons between A and Alist
%                                             % and if wished return them in struct <results>
%
%   ================ % initiate comparison with or without Alist
%
%   M = analyse.CompareModels(A); % Initiate with true network
%   % or
%   M = analyse.CompareModels(A,Alist); % Initiate and do comparison returning an object
%
%   ================ % query object
%
%   results = struct(M)       % Returns struct with measure-named fields and values
%   <measure> = M.<measure>;  % Returns array of <measure> values
%
%   results = compare(M,Alist,selected)  % returns the last results acquired with Alist as a struct;
%
%   ================ % get meaasures in a nice form
%
%   measures = properties(M);
%
%   ================ % find out what comparisons can be made
%
%   help analyse.CompareModels.show
%

    properties (Hidden = true)
        tol = eps;    % Link strength tolerance
    end

    properties (SetAccess = public, AbortSet = true)
        A           % Gold standard network
    end

    properties (Hidden = true, SetAccess = private)
        SA           % Singular values of true A
        UA           % Left singular vectors of true A
        VA           % Right singular vectors of true A
        LA           % Eigenvalues of true A
        QA           % Eigenvectors of true A
        DGA          % Adjacency matrix of A
        STA          % Signed topology of A
        N            % # Nodes
        ntl          % # True Links
        npl          % # possible links
    end

    properties (SetAccess = private)
        %% System measures

        abs2norm    % Absolute induced 2-norm
        rel2norm    % Relative induced 2-norm
        maee        % Max absolute element error
        mree        % Max relative element error
        mase        % Max absolute singular value error
        mrse        % Max relative singular value error
        masde       % Max absolute singular direction error
        mrsde       % Max relative singular direction error
        maeve       % Max absolute eigen value error
        mreve       % Max relative eigen value error
        maede       % Max absolute eigen direction error
        mrede       % Max relative eigen direction error
        afronorm    % Absolute Frobenius norm, equivalent to 2-norm of A vectorized
        rfronorm    % Relative Frobenius norm
        al1norm     % l1-norm of zero elements
        rl1norm     % Relative l1-norm of zero elements
        n0larger    % # zero elements larger than smallest nonzero element of A
        r0larger    % # zero elements larger than smallest nonzero element of A/# zero elements in A

        %% Signed Topology measures

        ncs         % # Correct signs
        sst         % Similarity of signed topology
        sst0        % Similarity of signed topology of nonzero elements of A

        %% Correlation measures

        plc         % Pearson's linear correlation coefficient

        %% Graph measures

        nlinks      % # Links in estimated network
        TP          % # True Positives
        TN          % # True Negatives
        FP          % # False Positives
        FN          % # False Negatives
        sen         % Sensitivity TP/(TP+FN)
        spe         % Specificity TN/(TN+FP)
        comspe      % Complementary specificity 1-Specificity
        pre         % Precision TP/(TP+FP)
        TPTN        % Number of links that is present and absent in both networks (TP+TN)
        structsim   % Structural similarity (TP+TN)/#Nodes^2
        MCC         % Matthews correlation coefficient

        %% Directed graph measures

        TR          % True Regulation
        TZ          % True Zero
        FI          % False Interaction
        FR          % False Regulation
        FZ          % False Zero
        dirsen      % Directed sensitivity
        dirspe      % Directed specificity
        dirprec     % Directed precision
        SMCC        % Signed Matthews correlation coefficient
    end

    methods
        function M = CompareModels(varargin)
            if nargin >= 1
                M.A = varargin{1};
            end

            if nargin >= 2
                Alist = varargin{2};
                M = compare(M,varargin{2:end});
            end
        end

        function M = set.A(M,net)
            A = M.A;
            if ~isempty(A) && ~isempty(net)
                warning('True A already set, overwriting')
            end
            if isa(net,'datastruct.Network')
                if length(size(net.A)) > 2
                    error('3d matrices are not allowed as golden standard')
                end
                M.A = net.A;
            else
                if length(size(net)) > 2
                    error('3d matrices are not allowed as golden standard')
                end
                M.A = net;
            end
            M = setA(M,M.A);
        end

        function M = setA(M,net)
        % Calculate network properties
            M.A = net;
            [UA SA VA] = svd(M.A);
            if issquare(M.A)
                [QA LA] = eig(M.A);
                LA = diag(LA);
                [crap index] = sort(abs(LA),1,'descend');
                LA = LA(index);
                QA = QA(:,index);
                M.LA = LA;
                M.QA = QA;
            end
            Z = zeros(size(M.A));
            DiGraphA = logical(Z);
            DiGraphA(abs(M.A) > M.tol) = true;
            STopoA = Z;
            STopoA(M.A > M.tol) = 1;
            STopoA(M.A < -M.tol) = -1;
            M.UA = UA;
            M.SA = SA;
            M.VA = VA;
            M.DGA = DiGraphA;
            M.STA = STopoA;
            M.N = size(M.A,1);
            M.npl = prod(size(M.A));
            M.ntl = sum(sum(DiGraphA));
        end

        function M = delA(M)
        % Deletes the compare network when object is no longer associated with a single network
            M = setA(M,[])
        end

        function zsst = sst2z(M)
        % Calculates similarity of signed topology for an empty network
            zsst = sum(sum(M.A == 0))/M.npl;
        end

        function M = system_measures(M,Alist)
        % Calculate only system measures
        % system_measures(M,Alist)

            for i=1:size(Alist,3)
                T = Alist(:,:,i);
                M.abs2norm(length(M.abs2norm)+1) = norm(T-M.A);
                try
                    M.rel2norm(length(M.rel2norm)+1) = norm(pinv(T)*T-M.A);
                catch err
                    M.rel2norm(length(M.rel2norm)+1) = NaN;
                end
                M.maee(length(M.maee)+1) = max(max(abs(T-M.A)));
                M.mree(length(M.mree)+1) = max(max(diag(1./max(abs(T),[],2))*abs(T-M.A)));

                S = svd(T);
                M.mase(length(M.mase)+1) = max(abs(S-diag(M.SA)));
                M.mrse(length(M.mrse)+1) = max(abs(diag(S)-M.SA)*(1./S));

                temp = M.UA'*T*M.VA;
                M.masde(length(M.masde)+1) = max(abs(diag(temp)-diag(M.SA)));
                M.mrsde(length(M.mrsde)+1) = max(diag(abs(diag(temp)-diag(M.SA)))*(1./diag(temp)));

                L = eig(T);
                [crap index] = sort(abs(L),1,'descend');
                L = L(index);
                M.maeve(length(M.maeve)+1) = max(abs(L-M.LA));
                M.mreve(length(M.mreve)+1) = max(abs(diag(1./L)*(L-M.LA)));

                temp = pinv(M.QA)*T*M.QA;
                temp2 = abs(diag(temp)-M.LA);
                M.maede(length(M.maede)+1) = max(temp2);
                M.mrede(length(M.mrede)+1) = max(abs(diag(temp2)*(1./diag(temp))));

                M.afronorm(length(M.afronorm)+1) = norm(T-M.A,'fro');
                try
                    M.rfronorm(length(M.rfronorm)+1) = norm(pinv(T)*(T-M.A),'fro');
                catch err
                    M.rfronorm(length(M.rfronorm)+1) = NaN;
                end
                M.al1norm(length(M.rl1norm)+1) = sum(abs(T(~M.DGA)));
                M.rl1norm(length(M.rl1norm)+1) = sum(abs(T(~M.DGA)))/sum(sum(abs(M.A)));
                M.n0larger(length(M.n0larger)+1) = sum(abs(T(~M.DGA)) > min(abs(M.A(M.DGA))));
                M.r0larger(length(M.r0larger)+1) = sum(abs(T(~M.DGA)) > min(abs(M.A(M.DGA))))/sum(sum(~M.DGA));
            end
        end

        function M = topology_measures(M,Alist)
        % Calculate only topololigcal measures
        % topology_measures(M,Alist)

            for i=1:size(Alist,3)
                T = Alist(:,:,i);
                STopoT = zeros(size(M.A));
                STopoT(T > M.tol) = 1;
                STopoT(T < -M.tol) = -1;
                temp = M.STA == STopoT;

                M.ncs(length(M.ncs)+1) = sum(sum(temp));
                M.sst(length(M.sst)+1) = sum(sum(temp))/M.npl;
                M.sst0(length(M.sst0)+1) = sum(sum(temp(M.DGA)))/M.ntl;
            end
        end

        function M = correlation_measures(M,Alist)
        % calculate only correlation measures
        % correlation_measures(M,Alist)

            for i=1:size(Alist,3)
                T = Alist(:,:,i);
                M.plc(length(M.plc)+1) = corr(M.A(:),T(:));
            end
        end

        function M = graph_measures(M,Alist)
        % Calculate only non directional graph measures
        % graph_measures(M,Alist)

            Z = zeros(size(M.A));
            for i=1:size(Alist,3)
                T = Alist(:,:,i);
                DiGraphT = logical(Z);
                DiGraphT(abs(T) > M.tol) = true;

                M.nlinks(length(M.nlinks)+1) = nnz(DiGraphT);
                M.TP(length(M.TP)+1) = sum(sum(and(M.DGA,DiGraphT)));
                M.TN(length(M.TN)+1) = sum(sum(and(~M.DGA,~DiGraphT)));
                M.FP(length(M.FP)+1) = sum(sum(and(~M.DGA,DiGraphT)));
                M.FN(length(M.FN)+1) = sum(sum(and(M.DGA,~DiGraphT)));
                M.sen(length(M.sen)+1) = M.TP(end)/(M.TP(end)+M.FN(end));
                M.spe(length(M.spe)+1) = M.TN(end)/(M.TN(end)+M.FP(end));
                M.comspe(length(M.comspe)+1) = M.FP(end)/(M.TN(end)+M.FP(end));
                M.pre(length(M.pre)+1) = M.TP(end)/(M.TP(end)+M.FP(end));
                M.TPTN(length(M.TPTN)+1) = M.TP(end) + M.TN(end);
                M.structsim(length(M.structsim)+1) = M.TPTN(end)/M.npl;

                n = (M.TP(end) + M.FP(end)) * (M.TP(end) + M.FN(end)) * (M.TN(end)+M.FP(end)) * (M.TN(end)+M.FN(end));
                if n == 0
                    M.MCC(length(M.MCC)+1) = 0;
                else
                    M.MCC(length(M.MCC)+1) = (M.TP(end)*M.TN(end)-M.FP(end)*M.FN(end))/sqrt(n);
                end
            end
        end

        function M = dirGraph_measures(M,Alist)
        % Calculate only directional graph measures
        % dirGraph_measures(M,Alist)

            for i=1:size(Alist,3)
                T = Alist(:,:,i);
                STT = zeros(size(M.A));
                STT(T > M.tol) = 1;
                STT(T < -M.tol) = -1;

                M.TR(length(M.TR)+1) = nnz( M.STA == 1 & STT == 1) + nnz(M.STA == -1 & STT == -1);
                M.TZ(length(M.TZ)+1) = nnz( and(~M.STA,~STT) );
                M.FI(length(M.FI)+1) = nnz( M.STA == 1 & STT == -1) + nnz(M.STA == -1 & STT == 1 );
                M.FR(length(M.FR)+1) = nnz( and(~M.STA,abs(STT)) );
                M.FZ(length(M.FZ)+1) = nnz( and(abs(M.STA),~STT) );

                if (M.TR(end) + M.FI(end) + M.FZ(end)) == 0
                    M.dirsen(length(M.dirsen)+1) = 0;
                else
                    M.dirsen(length(M.dirsen)+1) = M.TR(end)/(M.TR(end) + M.FI(end) + M.FZ(end));
                end

                if (M.TZ(end) + M.FR(end)) == 0
                    M.dirspe(length(M.dirspe)+1) = 0;
                else
                    M.dirspe(length(M.dirspe)+1) = M.TZ(end)/(M.TZ(end) + M.FR(end));
                end

                if (M.TR(end) + M.FI(end) + M.FR(end)) == 0
                    M.dirprec(length(M.dirprec)+1) = 0;
                else
                    M.dirprec(length(M.dirprec)+1) = M.TR(end)/(M.TR(end) + M.FI(end) + M.FR(end));
                end

                n = (M.TR(end) + M.FR(end)) * (M.TR(end) + M.FI(end) + M.FZ(end)) * (M.TZ(end) + M.FR(end)) * (M.TZ(end) + M.FI(end) + M.FZ(end));
                if n == 0
                    M.SMCC(length(M.SMCC)+1) = 0;
                else
                    M.SMCC(length(M.SMCC)+1) = (M.TR(end)*M.TZ(end) - M.FR(end)*(M.FI(end) + M.FZ(end))) / sqrt(n);
                end

            end
        end

        function varargout = compare(M,Alist,varargin)
        % Compares a 3D array of networks to the true network and returns the
        % selected measures
        %
        % <results => compare(obj, Alist [, selected,diagonal])
        %
        % Input variables:
        % ===============
        %
        % obj:      CompareModels object
        % Alist:    A list of networks to compare
        % selected: Only used if output is defined then it will return
        %           only the selected measures. Should be a vector of
        %           indices corresponding to indices given in show(obj).
        %
        % Output variables:
        % ================
        %
        %       M:  An object with all measures
        % results:  A struct with selected measures if any.

            A = M.A;
            if isempty(A)
                error('True network, A, needs to be set')
            end

            if nargin == 3
                selected = varargin{1};
            end

            % Eigenvalues can not be calculated for a non square matrix
            if issquare(M.A)
                M = system_measures(M,Alist);
                M = topology_measures(M,Alist);
                M = correlation_measures(M,Alist);
                M = graph_measures(M,Alist);
                M = dirGraph_measures(M,Alist);
            else
                M = topology_measures(M,tools.rmdiag(Alist));
                M = correlation_measures(M,tools.rmdiag(Alist));
                M = graph_measures(M,tools.rmdiag(Alist));
                M = dirGraph_measures(M,tools.rmdiag(Alist));
            end

            if nargout == 1 & nargin == 3
                allprops = show(M);
                results = struct([]);
                if exist('selected','var')
                    for i=selected
                        results(1).(allprops{i}) = M.(allprops{i});
                    end
                end
                varargout{1} = results;
            else
                varargout{1} = M;
            end
        end

        function varargout = show(M,varargin)
        % Get a list of measures or with input string get the index of that
        % measure if output argument is given, otherwise it will print
        % to terminal.
        %
        % {measures} = show(M [,name/number])
        %
            tmp = properties(M);
            allprops = tmp(find(~strcmp(tmp,'A')));
            % allprops = fieldnames(rmfield(tmp,'A'));

            if nargout == 0
                for i=1:length(allprops)
                    fprintf('\t%i\t%s\n',i,allprops{i})
                end
            elseif nargout == 1
                if length(varargin) == 0
                    varargout{1} = allprops;
                else
                    measure = varargin{1};
                    if isa(measure,'char')
                        varargout{1} = find(strcmp(allprops,measure));
                    elseif isa(measure,'double')
                        varargout{1} = allprops{measure};
                    end
                end
            end
        end

        function varargout = max(M,varargin);
        % Find maximum values over rows of specified measures or all values dependant on
        % a single measure.
        % [maximums, maxind] = max(M [,<measure>])
        %
        %
        %   Input Arguments: max(M [,<measure>])
        %   ================
        %   M         = Measure object
        %   <measure> = Optional measure to get maximum value of as given by the
        %               measures availible by show(M), may be both string or integer. ('all')
        %
        %   Output Arguments: [maximums, maxind]
        %   ================
        %   maximums  = If <measure> is not given then returns maximum values for each measure
        %               else the maximum value of given <measure> along with corresponding
        %               value for all other measures.
        %   maxind    = Gives the index of the maximum value for the given measures.
        %

            m = 'all';
            if nargin >= 1
                props = show(M);
            end
            if nargin == 2;
                m = varargin{1};
            end
            if isa(m,'double')
                m = props{m};
            end

            maximums = [];
            if strcmp(m,'all')
                for i=1:length(props)
                    [maximums(i,:), maxind(i,:)] = max((M.(props{i}))',[],1);
                end
            else
                ind = find(strcmp(props,m));
                [junk,maxind] = max((M.(props{ind}))',[],1);
                for i=1:length(props)
                    prop = M.(props{i})';
                    for j=1:size(prop,2)
                        maximums(i,j) = prop(maxind(j),j);
                    end
                end
            end

            varargout{1} = maximums;

            if nargout == 2
                varargout{2} = maxind;
            end
        end

        function varargout = maxmax(M,varargin);
        % Function returning a new CompareModels object with max for all measures.

            T = analyse.CompareModels();
            props = show(T);

            if length(varargin) == 0
                maxes = max(M);
                for i=1:length(props)
                    T.(props{i}) = maxes(i,:)';
                end
            else
                for j=1:length(varargin)
                    mess = varargin{j};
                    maxes = max(M,mess);
                    for i=1:length(props)
                        T.(props{i}) = [T.(props{i}), maxes(i,:)];
                    end
                end
            end
            varargout{1} = T;
        end

        function varargout = min(M,varargin);
        % Find minimum values over rows of specified measures or all values dependant on
        % a single measure.
        %
        %   Input Arguments: min(M [,<measure>])
        %   ================
        %   M         = Measure object
        %   <measure> = Optional measure to get minimum value of as given by the
        %               measures availible by show(M), may be both string or integer. ('all')
        %
        %   Output Arguments: [minimums, minind]
        %   ================
        %   minimums  = If <measure> is not given then returns minimum values for each measure
        %               else the minimum value of given <measure> along with corresponding
        %               value for all other measures.
        %   minind    = Gives the index of the minimum value for the given measures.
        %

            m = 'all';
            if nargin >= 1
                props = show(M);
            end
            if nargin == 2;
                m = varargin{1};
            end
            if isa(m,'double')
                m = props{m};
            end

            minimums = [];
            if strcmp(m,'all')
                for i=1:length(props)
                    [minimums(i,:), minind(i,:)] = min((M.(props{i}))',[],1);
                end
            else
                ind = find(strcmp(props,m));
                [junk,minind] = min((M.(props{ind}))',[],1);
                for i=1:length(props)
                    prop = M.(props{i})';
                    for j=1:size(prop,2)
                        minimums(i,j) = prop(minind(j),j);
                    end
                end
            end

            varargout{1} = minimums;

            if nargout == 2
                varargout{2} = minind;
            end
        end

        function varargout = minmin(M,varargin);
        % Function returning a new NetworkCompare object with min for all measures.
            T = analyse.CompareModels();
            props = show(T);

            if length(varargin) == 0
                mines = min(M);
                for i=1:length(props)
                    T.(props{i}) = mines(i,:)';
                end
            else
                for j=1:length(varargin)
                    mess = varargin{j};
                    mines = min(M,mess);
                    for i=1:length(props)
                        T.(props{i}) = [T.(props{i}), mines(i,:)];
                    end
                end
            end
            varargout{1} = T;
        end

        function pM = plus(M,N)
        % Execute plus operation for each measure between each measure.
            props = show(M);
            pM = analyse.CompareModels(M.A);
            for i=1:length(props)
                pM.(props{i}) = M.(props{i}) + N.(props{i});
            end
        end

        function M = vertcat(M,varargin)
        % Use vertical concatenation
            props = show(M);
            for j=1:length(varargin)
                N = varargin{j};
                for i=1:length(props)
                    M.(props{i}) = [M.(props{i}); N.(props{i})];
                end
            end
        end

        function M = horzcat(M,varargin)
        % Use horizontal concatenation
            props = show(M);
            for j=1:length(varargin)
                N = varargin{j};
                for i=1:length(props)
                    M.(props{i}) = [M.(props{i}), N.(props{i})];
                end
            end
        end

        function transM = transpose(M)
        % transposes all measures in object
            props = show(M);
            for i=1:length(props)
                M.(props{i}) = M.(props{i})';
            end
        end

        function M = stack(M,varargin)
        % stack each measure in the 3rd dimension (concatenation in 3rd dim)
            props = show(M);
            for j=1:length(varargin)
                N = varargin{j};
                for i=1:length(props)
                    M.(props{i}) = cat(3,M.(props{i}),N.(props{i}));
                end
            end
        end

        function mM = mean(M)
        % calculate mean of all measures,
        % if matrix then columnwise operation.
            props = show(M);
            mM = analyse.CompareModels(M.A);
            for i=1:length(props)
                mM.(props{i}) = mean(M.(props{i}));
            end
        end

        function vM = var(M)
        % calculate variance of all measures,
        % if matrix then columnwise operation.
            props = show(M);
            vM = analyse.CompareModels(M.A);
            for i=1:length(props)
                vM.(props{i}) = var(M.(props{i}));
            end
        end

        function T = getIndex(M,index)
        % returns an object with measure at index for each measure.
            T = analyse.CompareModels();
            props = show(T);

            for i=1:length(props)
                for j=1:size(M.(props{i}),1)
                    tmp(j) =  M.(props{i})(j,index(j));
                end
                T.(props{i}) = [T.(props{i}), tmp'];
            end
        end

        function T = reshape(M,varargin)
        % returns an object with reshaped each measure
            T = analyse.CompareModels();
            props = show(T);
            for i=1:length(props)
                T.(props{i}) = reshape(M.(props{i}),varargin{:});
            end
        end

        function square = issquare(M,A)
        % helper function if issquare
            square = true;
            [n,m] = size(A);

            if n ~= m
                square = false;
            end
        end
    end
end
