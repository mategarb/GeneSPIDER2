function [Y, E, stdE] = scnoise(X, P, dropout, zinbr, SNR)
% Create a noise matrix E based on log normal and zinb distributions
% with preferential attachment
%
% A = datastruct.scalefree(N,n[, seed])
%
% X:    noise free data
% P:    perturbetion matrix
% dropout:  droput rate 0-1, defining how many zeros in the data
% zinb: proportion of zinb noise 0-1, 0 means only log normal noise is used

s = svd(X);
stdE = s(size(P,1))/(SNR*sqrt(chi2inv(1-analyse.Data.alpha,numel(P))));
lnE = stdE*log2(random('LogNormal',0,1,size(P))); % noise model for lognormal
nbE = stdE*log2(random('Negative Binomial',20,0.5,size(P))/20); % negative binomial noise
ziE = random('Binomial',1,1-dropout,size(P)); % zero inflation, p is the proportion of 0s
zinbE = nbE.*ziE; %zinb noise

[rows, cols] = find(P~=0);
inds = [rows cols];

ngzinb = round(length(unique(inds(:,1)))*zinbr); % number of genes where zinb is applied

Y = X;
Y(:,inds(ismember(inds(:,1), 1:ngzinb)==1,2)) = (X(:,inds(ismember(inds(:,1), 1:ngzinb)==1,2)) + nbE(:,inds(ismember(inds(:,1), 1:ngzinb)==1,2))).*ziE(:,inds(ismember(inds(:,1), 1:ngzinb)==1,2));
Y(:,inds(ismember(inds(:,1), 1:ngzinb)==0,2)) = X(:,inds(ismember(inds(:,1), 1:ngzinb)==0,2)) + lnE(:,inds(ismember(inds(:,1), 1:ngzinb)==0,2));
E = zinbE;
E(:,inds(ismember(inds(:,1), 1:ngzinb)==0,2)) = lnE(:,inds(ismember(inds(:,1), 1:ngzinb)==0,2));

end