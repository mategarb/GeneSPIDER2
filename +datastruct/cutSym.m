function A = cutSym(A,inchance,feedback_fraction)
% removes links in a symmetric matrix with a probability to be in or out degree.
% A small chance of keeping feedback loops.
%
% A = datastruct.cutSym(A,inchance,feedback_fraction)
%
% A                 : network matrix nxn where i=outgoing, j=incoming
% inchance          : the chance of removing in i instead of j
% feedback_fraction : feedback = inchance/X where feedback_fraction=X. (default = 20)

N = size(A,1);

if ~exist('feedback_fraction','var')
    feedback = inchance/20;
else
    feedback = inchance/feedback_fraction;
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
