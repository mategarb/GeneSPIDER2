function A = scalefree(N,sparsity)
% Create a scalefree network with N nodes and specific sparseness
% with preferential attachment
%
% A = datastruct.scalefree(N,sparsity)
%
% A: undirected scalefree network matrix
% N: number of nodes
% sparsity: sparsity degree [0,1]

m0 = round(sparsity*N);
seed = logical(full(sprand(m0*2,m0*2-1,m0/m0^2)));
k = 0;
while rank(double(seed)) < min(size(seed))
    seed = randn(size(seed)).*logical(full(sprand(m0*2,m0*2-1,m0/m0^2)));
    % seed(floor(rand*numel(seed)+1)) = 1;
    k = k + 1;

    if ~mod(k,100)
        fprintf('k = %d\n',k)
    end
end

tmp = zeros(m0*2);
for i=1:size(seed,1)
    tmp(i,i+1:end) = seed(i,i:end);
end

for i=2:size(seed,1)
    tmp(i,1:i-1) = seed(i,1:i-1);
end
seed = logical(tmp+tmp');

A = zeros(N);
A(1:size(seed,1),1:size(seed,2)) = seed;

for i=(m0*2+1):N

    if rand < rem(sparsity,floor(sparsity))
        m = ceil(m0);
    else
        m = floor(m0);
    end
    if m == 0
        m = 1;
    end

    k = 0;
    while k < m
        ps = 0;
        r = rand;
        for inode = 1:i-1
            pl = sum(A(inode,1:i-1))/nnz(A(1:i-1,1:i-1));
            ps = ps + pl;
            if r < ps
                if A(i,inode) ~= 1
                    A(i,inode) = 1;
                    A(inode,i) = 1;
                    k = k + 1;
                    break
                else
                    break
                end
            end
        end
    end
end
