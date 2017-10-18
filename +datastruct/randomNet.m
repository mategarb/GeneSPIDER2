function A = randomNet(N,n)
% creates a random network with N nodes and specific sparseness with no self loops.
%
% A = datastruct.randomNet(N,n)
%
% N: number of nodes
% n: whished average number of links per node if n>=1, else relative sparsity of total possible links.
%
% A: undirected random network matrix

if n<1
    sparsity = n;
else
    sparsity = n/N;
end

A = zeros(N);
tspar = round(N^2*sparsity*N^2/(N*(N-1)))/N^2;
tmp = (1+rand*9)*full( sprandn(N,N-1,tspar) );
for i=1:size(tmp,1)
    A(i,i+1:end) = tmp(i,i:end);
end

for i=2:size(tmp,1)
    A(i,1:i-1) = tmp(i,1:i-1);
end
