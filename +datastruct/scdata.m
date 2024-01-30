function [Y, X, P, Ed, Eg] = scdata(A, P, options)
% Create a noise matrix E based on log normal and zinb distributions
% with preferential attachment

    % argument list with defaults
    arguments
        A double % GRN
        P double % perturbation matrix
        options.SNR (1,1) {mustBeNumeric} = 0.05; % signal-to-noise ratio
        options.raw_counts (1,1) {mustBeNumericOrLogical} = true % if true, raw counts are outputed, if false fold-change
        options.left_tail (1,1) {mustBeNumeric} = 1  % power of left tail in Pk distribution (>1 makes left tail longer)
        options.right_tail (1,1) {mustBeNumeric} = 2 % power of right tail in Pk distribution (>1 makes right tail longer)
        options.negbin_prob (1,1) {mustBeInRange(options.negbin_prob,0,1)} = 0.5 % mean for negative binomial used to simulate pseudo control
        options.disper (1,1) {mustBeNumeric} = 1 % dispersion parameter in Pk (Svensson et al.)
        options.n_clusts (1,1) {mustBeNumeric} = 5 % theoretical number of clusters
        options.logbase {mustBeNumeric} = 10 % log base to use in reversing FC, something that control perturbation strength
        options.ds_min {mustBeInRange(options.ds_min,0,1)} = 0.2 % minimum dissimilarity between clusters 
        options.ds_max {mustBeInRange(options.ds_max,0,1)} = 0.6 % maximum dissimilarity between clusters
    end  

%% prepare data + input noise
Net = datastruct.Network(A, 'scnet');
P = P(:,randperm(size(P,2))); % shuffling P in order to distribute expression across clusters
X = Net.G*P; %corresponding response Y, G is the static gain matrix (inverse of A (network matrix))

stdE = sqrt(var(X(:))/options.SNR);
Eg = stdE*randn(size(P)); % noise matrix
X = X + Eg;

%% generating theoretical clusters

% theoretical clusters
clst = options.n_clusts; % number of clusters
tots = size(P,2);
int1 = round(size(P,2)/clst);
invs = 0;
for it = 1:clst
    if tots > int1
        invs(it) = int1;
        tots = tots - int1;
    else
        invs(it) = tots;
    end
end

if sum(invs) < size(P,2)
    invs(end) = invs(end) + (size(P,2) - sum(invs));
end

Ab = {};
for z = 1:clst
    Ab{z} = ones(invs(z), 1);
end
bdg = blkdiag(Ab{:})';

% all possible theoretical clusters
npos = nchoosek(1:clst,2);
bdg2 = zeros(size(npos,1), size(P,2));
for u = 1:size(npos,1)
    bdg2(u,:) = sum(bdg(npos(u,:),:));
end
bdg3 = [bdg; bdg2];

%% generating pseudocontrol
psc = nbinrnd(1,options.negbin_prob,[1, size(P,1)]);  % pseduocontrol mean from negative binomial
psc = psc + rand(1, length(psc)); % add some 0-1 rand that will make it behave as non-discrete mean
psc = sort(psc,'ascend');
psc = psc.^options.right_tail; % to make the tail greater and leave 0s and 1s as 0s and 1s

X2 = ((options.logbase.^(X))); % reverse (pseudo)log
vX = var(X2,0,2); % take variance to relate it with mean
[varx, ~] = sort(vX);
[~, idx] = ismember(vX, varx);
psc2 = psc(idx); % mean expression sorted according to gene variance

sigma = sqrt(var(X(:))/options.SNR);

sC = zeros(size(P));
for i = 1:length(psc2)
        tcl = bdg3(randi([1,size(bdg3,1)],1),:); % select cluster
        tm = psc2(i);
        rr1 = rand; % providing some variation to the clusters
        rr2 = rand; % providing some variation to the clusters
        m1 = tm + tm*((options.ds_max-options.ds_min)*rr1 + options.ds_min); % generate two means distant by sf
        m2 = tm - tm*((options.ds_max-options.ds_min)*rr2 + options.ds_min);
        m2 = m2 + m2*(sum(tcl==0)/length(tcl)); % adjust mean for unequal length of clusters
        sC(i,tcl==1) = round(lognrnd(log(m1), log(sigma + 1), [1, sum(tcl==1)]));
        sC(i,tcl==0) = round(lognrnd(log(m2), log(sigma + 1), [1, sum(tcl==0)]));
end

X3 = round(X2.*sC); % multiply to get pseudo counts
%% inflate with dropouts

ziE = ones(size(P));
me = abs(mean(X3,2));
dspr = options.disper;
Pk0 = ((dspr.^-1)./(me + dspr.^-1)).^(dspr.^-1); % probability of dropouts where k = 0, from the paper
Pk0(Pk0 > 1) = 1; % because some probs are slighyly greater than 1 (rounding error most probably)

ziE = ones(size(P));
for ij = 1:length(Pk0)
     ziE(ij,:) = random('Binomial', 1, (1-Pk0(ij)), [1, size(P,2)]); % zero inflation, p is the proportion of 0s
end

%% final conversion to ouput
if options.raw_counts
    Y = X3.*ziE; % inflate with dropouts
    Y(Y<0) = 0; % get rid of negatives
else
    Y = X3.*ziE;
    Y = abs(Y./sC); % going back to fold change
    Y(ismissing(Y)) = 0; % NAs from 0 division
    Y(Y~=0) = Y(Y~=0) + eps; % add small value to preserve log10 links with 1 FC (from rounding error)
    Y = log(Y)/log(options.logbase); % going back to log 10 fold change
    Y(isinf(Y)) = 0; % infs from log10(1)
end

ziE3 = ziE;
ziE3(X3 == 0) = 0;
Ed = ziE3; % final dropouts

end