function A = randomNet(N,sparsity)
% creates a random network with N nodes and specific sparseness.
%
% A = GeneSpider.randomNet(N,sparsity)
%
% A: undirected random network matrix
% N: number of nodes
% sparsity: sparsity degree [0,1]

A = zeros(N);
tmp = (1+rand*9)*full( sprandn(N,N-1,round((sparsity*N^2/(N*(N-1)))*N^2)/N^2) );
for i=1:size(tmp,1)
    A(i,i+1:end) = tmp(i,i:end);
end

for i=2:size(tmp,1)
    A(i,1:i-1) = tmp(i,1:i-1);
end

A = GeneSpider.stabalize(A);
for i=1:n
    for j=1:n
        if abs(A(i,j)) < tol || isnan(A(i,j))
            A(i,j) = 0;
        end
    end
end
