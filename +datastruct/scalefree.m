function A = scalefree(N,n,varargin)
% Create a scalefree network with N nodes and specific sparseness
% with preferential attachment
%
% A = datastruct.scalefree(N,n[, seed])
%
% N:    number of nodes
% n:    whished average number of links per node if n>=1, else relative sparsity of total possible links.
% seed: matrix to work as a seed, assumed to have size < N
%
% A:    undirected scalefree network matrix

rank_check = true;
if length(varargin) > 0
    for i=1:length(varargin)
        if isa(varargin{i},'logical')
            rank_check = varargin{i};
        else
            seed = varargin{i};
        end
    end
end

if n<1
    sparsity = n;
else
    sparsity = n/N;
end

m0 = round(sparsity*N);

if ~exist('seed','var')
    seed = logical(full(sprand(m0*2,m0*2-1,m0/m0^2)));
    k = 0;
    if rank_check
        while rank(double(seed)) < min(size(seed))
            seed = randn(size(seed)).*logical(full(sprand(m0*2,m0*2-1,m0/m0^2)));
            % seed(floor(rand*numel(seed)+1)) = 1;
            k = k + 1;

            if ~mod(k,100)
                fprintf('k = %d\n',k)
            end
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
end

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
                A(i,inode) = 1;
                A(inode,i) = 1;
                k = k + 1;
                break
            end
        end
    end
end
