function A = cutSym(A,pin,pout)
% removes links in a symmetric matrix with a probability to be in or out degree.
% A small chance of keeping feedback loops.
%
% A = datastruct.cutSym(A,pin,pout)
%
% A                 : network matrix nxn where i=outgoing, j=incoming
% pin               : the chance of removing in i instead of j
% pout              : the chance of removing in j instead of i
%
% the probability of a feedback loop is then
% pfeedback         : 1 - (pin+pout)

N = size(A,1);

for i=1:N-1
    for j=i+1:N
        if A(i,j) == 1
            r = rand;
            if r <= pin
                A(i,j) = 0;
            elseif r <= pin+pout
                A(j,i) = 0;
            end
        end
    end
end
