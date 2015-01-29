function Mhat = standardize(M)
%
% Standard normalize the given matrix so that:
%       1. The sum of all values in each row = 0
%       2. The sum of the square of all values in each row = #columns
%
% In reality, for each row, subtract the mean and then divide by the
% standard deviation.
%
% M: raw matrix (n x m)
% 
% Mhat: standardized matrix (n x m)
%
m = size(M, 2);

mu = mean(M, 2);
sigma = std(M, 1, 2); % population standard deviation
Mu = repmat(mu, 1, m);
Sigma = repmat(sigma, 1, m);

Mhat = (M - Mu) ./ Sigma;