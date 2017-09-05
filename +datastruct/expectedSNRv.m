function SNRv = expectedSNRv(data)
% Calculation of the expected $\SNR_{\phi\, \normall}$
% SNRv = expectedSNRv(data)
%
%

N = data.N; % Number of genes in the network that SNRv is estimated for
M = data.M; % Number of experiments in the data set that SNRv is estimated for
P = data.P % True strength of perturbations is unknown, so each row of A will be relative to the strength of the perturbation
Y = response(data);

% Prep = data.Prep;
% Yrep = data.Yrep;

% RegressorLengths = normv(Y','vec');
ResponseLengths = normv(Y,'vec');

RegressorLength = (RegressorLengths)/sqrt(size(Y',1)); % The expected norm of the response in a gene in an experiment.

lambda = data.lambda(1:length(data.lambda)/2);

SNRv = RegressorLength*sqrt(M)./sqrt(lambda*chi2inv(1-data.alpha,M));
