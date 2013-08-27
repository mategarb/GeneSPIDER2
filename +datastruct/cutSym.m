function A = cutSym(A,inchance)
% removes links in a symmetric matrix with a probability to be in or out degree
%
% A = GeneSpider.cutSym(A,inchance)
%
% A: network matrix nxm n=outdegree, m=indegree
% inchance: the chanse of removing in n instead of m

N = size(A,1);

for i=1:N-1
    for j=i+1:N
        if A(i,j) == 1
            if rand < inchance
                A(i,j) = 0;
            else
                A(j,i) = 0;
            end
        end
    end
end
