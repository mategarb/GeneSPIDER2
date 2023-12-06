function A = cutSym(A,pupper,plower)
% removes links in a symmetric matrix with settable probabilities for upper and lower triangular parts.
%
% A = datastruct.cutSym(A,pupper,plower)
%
% A                 : network matrix nxn
% pupper            : the chance of removing i,j
% plower            : the chance of removing j,i
%

N = size(A,1);

for i=1:N-1
    for j=i+1:N
        if A(i,j) ~= 0 && A(j,i) ~= 0
            r = rand;
            if r <= pupper
                A(i,j) = 0;
            elseif r <= pupper+plower
                A(j,i) = 0;
            end
        end
    end
end
