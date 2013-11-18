function candidateP = weightP(A,P)

[UA SA VA] = svd(A);
G = -pinv(A);
level = [0.1 0.01 0.01 0.001]; % Original. Why two of the same
tol = 10^-4; % Tolerance
limitSVY = 0.33; % Limit for how much the SVs of Y may differ from 1 (wished value)


k = 1;
removedelement = false;
while k <= length(level)

    %% Tries to optimise non-zero elements of P
    % P(abs(P) < level(end)) = 0;
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
                    % SVYsumplus = norm(SVY-ones(size(SVY)),1); % Sum of SVs difference from one
                    SVYsumplus = norm(log(SVY)-ones(size(SVY)),2); % Sum of SVs difference from one, this is not actually the sum

                    Ptemp(i,j) = Ptemp(i,j) - 2*level(k);
                    SVY = svd(G*Ptemp);
                    % SVYsumminus = norm(SVY-ones(size(SVY)),1); % Sum of SVs difference from one
                    SVYsumminus = norm(log(SVY)-ones(size(SVY)),2); % Sum of SVs difference from one, this is not actually the sum

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

    k = k+1; % Advance to next level
end
%% Store the resulting P
SVY = svd(G*P);
SVYeffall = max(abs(SVY-ones(size(SVY)))); % Largest difference in singular value from one
P(abs(P) < tol) = 0;
P = round(P*1/tol)*tol;
candidateP = P;
fprintf([datestr(clock) '\tSVYeff:' num2str(SVYeffall),'\n'])
