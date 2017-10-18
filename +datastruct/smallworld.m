function [A,s] = smallworld(N,k,p,undirected,varargin)
% generate a small world network
%
% [A,s] = smallworld(N,k,p,undirected,d)
%
% A: undirected or directed adjacency matrix
% s: random number generator state
%
% N:          number of nodes
% k:          inital connections for each node, assumed even
% p:          probability of rewiring
% undirected: logical if undirected or directed
% d:          probability of rewiring direction
%

if isempty(varargin)
    d = 0.5;
else
    d = varargin{1};
end

A = zeros(N,N);

s = rng();

for i=1:N
    for j=i+1:i+k/2
        t = j;
        if t > N
            t = t - N;
        end
        A(i,t) = 1;
        A(t,i) = 1;
        r = rand();
        if p > r
            notselection = mod([i-k/2:i+k/2],N)+1;
            selection = [1:N];
            selection(notselection) = 0;
            selection = selection(find(selection));
            l = randi(length(selection),1);
            sel = selection(l);
            while sel == t
                l = randi(length(selection),1);
                sel = selection(l);
            end

            A(i,t) = 0;
            A(t,i) = 0;

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
                    notselection = mod([i-k/2:i+k/2],N)+1;
                    selection = [1:N];
                    selection(notselection) = 0;
                    selection = selection(find(selection));
                    h = randi(length(selection),1);
                    sel = selection(h);
                    A(sel,i) = 1;
                end
            end
        end
    end
end
