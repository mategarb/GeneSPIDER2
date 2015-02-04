function A = cutSym(A,inchance)
% removes links in a symmetric matrix with a probability to be in or out degree.
% A small chance of keeping feedback loops.
%
% A = datastruct.cutSym(A,inchance)
%
% A: network matrix nxm n=outdegree, m=indegree
% inchance: the chance of removing in n instead of m

N = size(A,1);

if inchance < 0.5
    feedback = inchance/20;
else
    feedback = (1-inchance)/20;
end

inchance = inchance - feedback;

for i=1:N-1
    for j=i+1:N
        if A(i,j) == 1
            r = rand;
            if r < inchance
                A(i,j) = 0;
            elseif r > inchance+feedback
                A(j,i) = 0;
            end
        end
    end
end
