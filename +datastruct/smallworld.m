function [A,s] = smallworld(N,k,p,undirected,varargin)
% generate a small world network
%
% [A,s] = smallworld(N,k,p,undirected,d)
% 
% A: undirected or directed adjacency matrix
% s: random number generator state
% 
% N: number of nodes
% k: inital connections for each node, assumed even
% p: probability of rewireing
% undirected: logical if undirected or directed
% d: probability of rewiring direction
%

if isempty(varargin)
    d = 0.5
else
    d = varargin{1};
end

A = zeros(N,N);

s = rng();

for i=1:N
    for j=mod(i+1,N):(mod(i+k/2,N)+1)

        A(i,j) = 1;
        A(j,i) = 1;

        
        r = rand();
        if p > r
            l = randi(N,1);
            while any(l == mod([i-k/2:i+k/2],N))
                if l == j, continue; end
                l = randi(N,1);
            end
            
            A(i,j) = 0;
            A(j,i) = 0;
            
            if undirected
                A(i,l) = 1;
                A(l,i) = 1;
            else
                r = rand();
                if d > r
                    A(i,l) = 1;
                    A(l,i) = 1;
                else
                    A(i,l) = 1;
                    h = randi(N,1);
                    while any(h == mod([i-k/2:i+k/2],N))
                        h = randi(N,1);
                    end
                    A(h,i) = 1;
                end
            end
        end
    end
end

