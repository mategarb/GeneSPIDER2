function candidateP = optimalPstructure(A,nm,sp)
% optimalP tries to find optimal P via convex optimization.
% 
% candidateP = optimalP(A,nm,sp)
% 
% A: Network
% nm: dimensions of P, size(P)
% 
% 
% 

n = nm(1);
m = nm(2);
G = -pinv(A);
[UA SA VA] = svd(A);
G = -pinv(A); % Static gain of system
R = GramSchmidtOrth(-1+2*rand(max([n m]))); % Random R as a start
R = R(:,1:n)';
P = UA*SA*R; % P based on economy size svd
P = round(1000.*P)./1000;
Y = G*P;

candidateP = [];
for i=1:m
    fit = glmnet(G,Y(:,i));
    beta = fit.beta;    
    nPert = sum(logical(beta));

    j = max(find(nPert == sp));
    if isempty(j)
        j = max(find(nPert == sp-1));
    end
    if isempty(j)
        j = max(find(nPert == sp+2));
    end
    
    jPert = find(beta(:,j));
    
    candidateP(:,i) = beta(:,j);
    % candidateP(jPert,i) = beta(jPert,end);
    
end




