function candidateP = optimalRandomP(A,nCandidates,nm)
% candidateP = optimalRandomP(A,nCandidates,size(P))
%
% optimalRandomP calculates an optimal perturbation design of desired
% size for a network response to be informative when influenced by noise.
%
% A: Network matrix
% nCandidates: How many candidates should be generated
% nm: Vector [n,m]. Where m is # nodes and m is # experiments.

n = nm(1);
m = nm(2);
[UA SA VA] = svd(A);
G = -pinv(A); % Static gain of system
limitSVY = 0.33; % Limit for how much the SVs of Y may differ from 1 (wished value)

% level = logspace(-4,-1,4);
% level = [0.1 0.01 0.001];
level = [0.1 0.01 0.01 0.001]; % Original. Why two of the same
tol = 10^-4; % Tolerance
keepOneElementOnEachColumn = true; % Used to avoid experiments where nothing is perturbed


for r = 1:nCandidates

    %% Generates new random perturbation matrix
    R = GramSchmidtOrth(-1+2*rand(max([n m]))); % Random R as a start
    R = R(:,1:n)';
    P = UA*SA*R; % P based on economy size svd
    P = round(1000.*P)./1000;
    nrpert(r) = sum(sum(P ~= 0)); % Number of non-zero elements in P

    %% Optimises the P to be as sparse as possible
    k = 1;
    removedelement = false;
    while k <= length(level)
        P(abs(P) < 0.001) = 0;
        SVY = svd(G*P);
        SVYsumP = norm(SVY-ones(size(SVY)),2); % Sum of SVs difference from one
        fprintf([datestr(clock) '\t' num2str(r) '\t' num2str(k) '\tnrpert:' num2str(nrpert(r)) '\tSVYsumP:' num2str(SVYsumP,10),'\n'])


        %% Tries to set as many elements to zero as possible
        while true,
            %% Find the effect on SVs of Y of setting each element one-by-one of P to zero
            for i = 1:size(P,1),
                for j = 1:size(P,2),
                    if P(i,j) ~= 0
                        Ptemp = P;
                        Ptemp(i,j) = 0;
                        SVY = svd(G*Ptemp);
                        SVYeff(i,j) = max(abs(SVY-ones(size(SVY)))); % Find SV of Y that differs most from 1 (wished value)
                    end
                end
            end
            SVYeff(P == 0) = Inf; % Excludes all elements that already are zero
            if keepOneElementOnEachColumn,
                SVYeff(:,sum(P == 0,1) == n-1) = Inf; % Excludes all columns that only have one non-zero element
            end
            %% Find index of element of P to set to zero with smallest effect on SV of Y
            [val ind] = min(reshape(SVYeff,numel(SVYeff),1));
            if val <= limitSVY,
                P(ind) = 0;
                removedelement = true;
            else
                break
            end
        end

        %% Tries to optimise non-zero elements of P
        P(abs(P) < level(end)) = 0;
        SVY = svd(G*P);
        SVYsumPold = norm(SVY-ones(size(SVY)),2); % Sum of SVs difference from one
        SVYsumP = SVYsumPold;
        while true % Go through all elements one by one
            for i = 1:size(P,1)
                for j = 1:size(P,2)
                    if P(i,j) ~= 0
                        % Make the element unchanged, plus or minus one level,
                        % depending on which gives smallest criteria
                        Ptemp = P;

                        Ptemp(i,j) = Ptemp(i,j) + level(k);
                        SVY = svd(G*Ptemp);
                        SVYsumplus = norm(SVY-ones(size(SVY)),2); % Sum of SVs difference from one

                        Ptemp(i,j) = Ptemp(i,j) - 2*level(k);
                        SVY = svd(G*Ptemp);
                        SVYsumminus = norm(SVY-ones(size(SVY)),2); % Sum of SVs difference from one

                        [val ind] = min([SVYsumP SVYsumminus SVYsumplus]);
                        SVYsumP = val;
                        switch ind,
                            case 1,
                              % Keep same P
                            case 2,
                                P(i,j) = P(i,j) - level(k);
                            case 3,
                                P(i,j) = P(i,j) + level(k);
                        end
                    end
                end
            end
            if SVYsumP - SVYsumPold < -tol,
                % Improvement took place continue looping
                SVYsumPold = SVYsumP;
            else
                % No improvement so break while loop
                break
            end
        end

        %% Determines if should move to next level or not
        if removedelement,
            nrpert(r) = sum(sum(P ~= 0));
            removedelement = false; % Stay on same level and optimise elements again
        else
            k = k+1; % Advance to next level
        end

        %% Store the resulting P
        nrpert(r) = sum(sum(P ~= 0)); % Number of non-zero elements in P
        SVY = svd(G*P);
        SVYeffall(r) = max(abs(SVY-ones(size(SVY)))); % Largest difference in singular value from one
        Pall(:,:,r) = P;
        fprintf([datestr(clock) '\t' num2str(r) '\t\tnrpert:' num2str(nrpert(r)) '\tSVYeff:' num2str(SVYeffall(r)),'\n'])

    end
end

Pall(abs(Pall) < tol) = 0;

candidateP = Pall;
