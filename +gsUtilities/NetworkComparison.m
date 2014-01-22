classdef NetworkComparison < hgsetget
% CompareNetworks calculates difference measures between weighted network
% adjacency matrices.
% For a nice overview do:
% doc tools.NetworkCopmarison
%
%   Input Arguments: NetworkComparison(A, Alist, selected)
%   ================
%   (none)     Initiate object for later use.
%   A        = The 'true' weighted adjacency matrix of the network, to which all
%              other network matrices are compared.
%   Alist    = A network matrix or 3d array of network matrices that are
%              compared to the 'true' one.
%
%   Output Arguments:
%   ================
%   obj = NetworkComparision object.
%
%   Use Examples:
%   ================
%
%   M = NetworkComparison();                  % initiate comparision
%   % set(M,'tol',value)                      % set tolerance {eps}
%   M.A = A;                                  % set true A
%   <results> = compare(M,Alist[,selected]);  % Update M with comparisons between A and Alist
%                                             % and if wished return them in struct <results>
%
%   ================ % initiate comparison with or without Alist
%
%   M = NetworkComparison(A); % Initiate with true network
%   % or
%   M = NetworkComparison(A,Alist); % Initiate and do comparison returning an object
%
%   ================ % querie object
%
%   results = get(M)          % Returns struct with measure-named fields and values
%   <measure> = M.<measure>;  % Returns array of <measure> values
%
%   results = compare(M,Alist,selected)  % retruns the last results aquired with Alist as a strcut;
%
%   ================ % get meaasures in a nice form
%
%   measures = struct2cell(get(M));
%   % or
%   measures = struct2array(compare(M,Alist,selected));
%
%   ================ % find out what comparisons can be made
%
%   help tools.NetworkComparison.show
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
        FEL        % Fraction of Existing Links (for RNI = FPEL)

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
        function M = NetworkComparison(varargin)
            if nargin >= 1
                setA(M,varargin{1});
            end

            if nargin >= 2
                Alist = varargin{2};
                compare(M,varargin{2:end});
            end
        end

        function set.A(M,net)
            A = M.A;
            if ~isempty(A) && ~isempty(net)
                warning('True A already set, overwriting')
            end
            if isa(net,'GeneSpider.Network')
                M.A = net.A;
            else
                M.A = net;
            end
        end

        function setA(M,net)
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

        function delA(M)
        % Deletes the compare network when object is no longer associated with a single network
            setA(M,[])
        end
        
        function zsst = sst2z(M)
        % Calculates similarity of signed topology for an empty network
            zsst = sum(sum(M.A == 0))/M.npl;
        end

        function systemMeasures(M,Alist)
            for i=1:length(Alist(1,1,:))
                T = Alist(:,:,i);
                M.abs2norm(length(M.abs2norm)+1) = norm(T-M.A);
                M.rel2norm(length(M.rel2norm)+1) = norm(pinv(T)*T-M.A);
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
                M.rfronorm(length(M.rfronorm)+1) = norm(pinv(T)*(T-M.A),'fro');
                M.al1norm(length(M.rl1norm)+1) = sum(abs(T(~M.DGA)));
                M.rl1norm(length(M.rl1norm)+1) = sum(abs(T(~M.DGA)))/sum(sum(abs(M.A)));
                M.n0larger(length(M.n0larger)+1) = sum(abs(T(~M.DGA)) > min(abs(M.A(M.DGA))));
                M.r0larger(length(M.r0larger)+1) = sum(abs(T(~M.DGA)) > min(abs(M.A(M.DGA))))/sum(sum(~M.DGA));
            end
        end

        function topologyMeasures(M,Alist)
            for i=1:length(Alist(1,1,:))
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

        function correlationMeasures(M,Alist)
            for i=1:length(Alist(1,1,:))
                T = Alist(:,:,i);
                M.plc(length(M.plc)+1) = corr(M.A(:),T(:));
            end
        end

        function graphMeasures(M,Alist)
            Z = zeros(size(M.A));
            for i=1:length(Alist(1,1,:))
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
                M.FEL(length(M.FEL)+1) = nnz(M.DGA & DiGraphT)/M.ntl;
                
                n = (M.TP(end) + M.FP(end)) * (M.TP(end) + M.FN(end)) * (M.TN(end)+M.FP(end)) * (M.TN(end)+M.FN(end));
                if n == 0
                    M.MCC(length(M.MCC)+1) = 0;
                else
                    M.MCC(length(M.MCC)+1) = (M.TP(end)*M.TN(end)-M.FP(end)*M.FN(end))/sqrt(n);
                end
            end
        end

        function dirGraphMeasures(M,Alist)
            for i=1:length(Alist(1,1,:))
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
        % obj:      NetworkComparison object
        % Alist:    A list of networks to compare
        % selected: Only used if output is defined then it will return
        %           only the selected measures. Should be a vector of
        %           indices corresponding to indices given in show(obj).
        %
        % Output variables:
        % ================
        %
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
                systemMeasures(M,Alist);
                topologyMeasures(M,Alist);
                correlationMeasures(M,Alist);
                graphMeasures(M,Alist);
                dirGraphMeasures(M,Alist);
            else
                topologyMeasures(M,tools.rmdiag(Alist));
                correlationMeasures(M,tools.rmdiag(Alist));
                graphMeasures(M,tools.rmdiag(Alist));
                dirGraphMeasures(M,tools.rmdiag(Alist));
            end

            if nargout == 1
                allprops = show(M);
                results = struct([]);
                if exist('selected','var')
                    for i=selected
                        results(1).(allprops{i}) = M.(allprops{i});
                    end
                end
                varargout{1} = results;
            end
        end

        function varargout = show(M,varargin)
        % Get a list of measures or with input string get the index of that
        % measure if output argument is given, otherwise it will print
        % to terminal.
        %
        % {measures} = show(M [,name/number])
        %
            tmp = set(M);
            allprops = fieldnames(rmfield(tmp,'A'));

            if nargout == 0
                for i=1:length(allprops)
                    fprintf('%i %s\n',i,allprops{i})
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
        % Find maximum values of specified measures or all values dependant on
        % a single measure.
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
                    [maximums(i,:), maxind] = max((M.(props{i}))');
                end
            else
                ind = find(strcmp(props,m));
                [junk,maxind] = max((M.(props{ind}))');
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
        % Function returning a new Compare object with max for all measures.
        end

        function varargout = min(M,varargin);
        % Find minimum values of specified measures or all values dependant on a single measure.
        %
        %   Input Arguments: min(M [,<measure>])
        %   ================
        %   M         = Measure object
        %   <measure> = Optional measure to get minimum value of as given by the
        %               measures availible by show(M), may be both string or integer. ('all')
        %
        %   Output Arguments: [minimums,minind]
        %   ================
        %   minimums  = If <measure> is not given then returns minimum values for each measure
        %               else the minimum value of given <measure> along with corresponding
        %               value for all other measures.
        %   maxind    = Gives the index of the minimum value for the given measures.
        %

            m = 'all';
            if nargin == 1
                props = show(M);
            end
            if nargin == 2;
                m = varargin{1};
            end
            if isa(m,'double')
                mainprop = props{m};
            end

            minimums = [];
            if strcmp(m,'all')
                for i=1:length(props)
                    [minimums(i,:),minind] = min(M.(props{i}));
                end
            else
                ind = find(strcmp(props,m));
                minind = find(M.(props{ind}) == min(M.(props{ind})));
                for i=1:length(props)
                    minimums(i,:) = M.(props{i})(minind);
                end
            end

            varargout{1} = minimums;

            if nargout == 2
                varargout{2} = minind;
            end
        end

        function pM = plus(M,N)
            props = show(M);
            pM = tools.NetworkComparison(M.A);
            for i=1:length(props)
                pM.(props{i}) = M.(props{i}) + N.(props{i});
            end
        end

        function vertM = vertcat(M,varargin)
            props = show(M);
            vertM = tools.NetworkComparison(M.A);
            for j=1:length(varargin)
                N = varargin{j};
                for i=1:length(props)
                    vertM.(props{i}) = [M.(props{i}); N.(props{i})];
                end
            end
        end

        function horzM = horzcat(M,varargin)
            props = show(M);
            horzM = tools.NetworkComparison(M.A);
            for j=1:length(varargin)
                N = varargin{j};
                for i=1:length(props)
                    propiM = M.(props{i});
                    propiN = N.(props{i});
                    if isrow(propiM), propiM = propiM'; end
                    if isrow(propiN), propiN = propiN'; end
                    horzM.(props{i}) = [propiM, propiN];
                end
            end
        end

        function mM = mean(M)
        % calculate mean of all measures
            props = show(M);
            mM = tools.NetworkComparison(M.A);
            for i=1:length(props)
                mM.(props{i}) = mean(M.(props{i}));
            end
        end

        function vM = var(M)
        % calculate variance of all measures
            props = show(M);
            mM = tools.NetworkComparison(M.A);
            for i=1:length(props)
                mM.(props{i}) = var(M.(props{i}));
            end
        end
    end
end
