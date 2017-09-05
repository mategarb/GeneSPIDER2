function [candidateP,Y] = optimalPstructure(A,nm,np)
% optimalP tries to find optimal P via convex optimization.
%
% candidateP = optimalP(A,nm,sp)
%
% A: Network
% nm: dimensions of P, size(P)
% np: number of  whiched perturbations / sample.
%

n = nm(1);
m = nm(2);
G = -pinv(A);
[UA SA VA] = svd(A);
G = -pinv(A); % Static gain of system
R = GramSchmidtOrth(-1+2*rand(max([n m]))); % Random R as a start
R = R(:,1:n)';
P = UA*SA*R; % P based on economy size svd
P = round(10000.*P)./1000;
Y = G*P;

candidateP = [];
for i=1:m
    sparse = np;
    fit = glmnet(G,Y(:,i));
    beta = fit.beta;
    nPert = sum(logical(beta));

    j = [];
    j = max(find(nPert == sparse));
    if isempty(j)
        j = max(find(nPert == sparse-1));
    end
    while isempty(j)
        sparse = sparse+1;
        j = max(find(nPert == sparse));
    end

    jPert = find(beta(:,j));

    candidateP(:,i) = beta(:,j);
    % candidateP(jPert,i) = beta(jPert,end);

end

nzr = find(sum(logical(candidateP)')==0);
for i=1:length(nzr)
    npert = sum(logical(candidateP));
    whatColumn = find(npert == min(npert));
    candidateP(nzr(i),whatColumn(1)) = -10+20*rand;
end
nzr = find(sum(logical(candidateP)')==0);
if ~isempty(nzr)
    warning('WTF')
end

rP = rank(candidateP);
if rP < min(m,n)
    warning('low rank candidate P')
    return
end

howMany = nnz(sum(logical(candidateP)));
if howMany ~= m
    candidateP = zeros(n,m);
end
