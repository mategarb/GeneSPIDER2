function varargout = DKL(p,q)
% Calculate the Kullback–Leibler divergence given two descrete probability distributions
% [kldiv, kla] = analyse.DKL(p,q)
%
% kldiv: is the Kullback–Leibler divergence
% kla: is the Kullback–Leibler Area that needs to be integrated.
% p: distribution of "true"
% q: distribution of "random"
% p and q must be of the same length and in the same order


% check if the probability p and q are probability distributions, i.e. sum to 1
if (abs(sum(p)) == 0) || (abs(sum(q)) == 0)
    error('p or q are not a probability distribution or contain no values.')
end

if (abs(sum(p) - 1) > 1e-10)
    p = p/sum(p);
end

if (abs(sum(q) - 1) > 1e-10)
    q = q/sum(q);
end

kla = p .* log2(p./q);
kldiv = sum(kla);


varargout{1} = kldiv;
if nargout > 1
    varargout{2} = kla;
end
