function [Y, E] = scnoise(X, P, options)
% Create a noise matrix E based on log normal and zinb distributions
% with preferential attachment

    % argument list with defaults
    arguments
        X double % noise-free data matrix
        P double % perturbation matrix
        options.SNR (1,1) {mustBeNumeric} = 0.5 % signal-to-noise ratio
        options.noise_model = "ZINB" % noise model, ZINB or ZIG (zero-inflated negative binomial or gaussian)
        options.raw_counts (1,1) {mustBeNumericOrLogical} = true % if true, raw counts are outputed, if false fold-change
        options.left_tail (1,1) {mustBeNumeric} = 1  % power of left tail in Pk distribution (>1 makes left tail longer)
        options.right_tail (1,1) {mustBeNumeric} = 1.5 % power of right tail in Pk distribution (>1 makes right tail longer)
        options.negbin_mean (1,1) {mustBeNumeric} = 0.5 % mean for negative binomial used to simulate pseudo control
        options.disper (1,1) {mustBeNumeric} = 1.5 % dispersion parameter in Pk (Svensson et al.)
        options.n_clusts (1,1) {mustBeNumeric} = 5 % theoretical number of clusters
        options.prob_clusts {mustBeInRange(options.prob_clusts,0,1)} = 0.6 % probability of gene belonging to a cluster
        options.prob_dev {mustBeInRange(options.prob_dev,0,1)} = 0.9 % probability of gene being deviated from Pk=0  
    end  

psc = nbinrnd(1,options.negbin_mean,[1, size(P,1)]); % pseduocontrol mean from negative binomial
psc = psc + rand(1, length(psc)); % add some 0-1 rand that will make it behave as non-discrete mean
psc = sort(psc,'ascend');
psc = psc.^options.right_tail; % to make the tail greater and leave 0s and 1s as 0s and 1s

if strcmp(options.noise_model, "ZIG") 

    s = svd(X);
    stdE = s(N)/(options.SNR*sqrt(chi2inv(1-analyse.Data.alpha,numel(P))));
    E = stdE*randn(size(P)); % noise matrix
    X = X + E;
end

X = ((10.^(X))); % reverse (pseudo)log
vX = var(X,0,2); % take variance to relate it with mean
[varx, ~] = sort(vX);
[~, idx] = ismember(vX, varx);
psc2 = psc(idx); % mean expression sorted according to gene variance

sC = zeros(size(P));
for i = 1:length(psc2)
        sC(i,:) = round(normrnd(psc2(i), 1, [1, size(P,2)])); % some pseudo control values to multiply
end

X = round(X.*sC); % multiply to get pseudo counts

ziE = ones(size(P));
me = abs(mean(X,2));
dspr = options.disper;
Pk0 = ((dspr.^-1)./(me + dspr.^-1)).^(dspr.^-1); % probability of dropouts where k = 0, from the paper
Pk0(Pk0 > 1) = 1; % because some probs are slighyly greater than 1 (rounding error most probably)

for i = 1:length(Pk0)
    ziE(i,:) = random('Binomial', 1, (1-Pk0(i))^options.left_tail, [1, size(P,2)]); % zero inflation, p is the proportion of 0s
end

% simulating UMAP clusters for random genes
clst = options.n_clusts;

Ab = {};
for z = 1:clst
    Ab{z} = ones(round(size(ziE,2)/clst), 1);
end
bdg = blkdiag(Ab{:})';

% all possible theoretical clusters
npos = nchoosek(1:clst,2);
bdg2 = zeros(size(npos,1), size(ziE,2));
for u = 1:size(npos,1)
    bdg2(u,:) = sum(bdg(npos(u,:),:));
end


ziE2 = ziE;
for xi = 1:size(ziE,1)
    rac = rand;
    if rac < options.prob_clusts
     tcl = bdg2(randi([1,size(bdg2,1)],1),:);
        if (sum(ziE(xi,:)) >= sum(tcl)) % more values than in theoretical
            new_tcl = tcl;
            ind_zer = find(new_tcl == 0);
            remo = round((sum(ziE(xi,:)) - sum(tcl))); % number of remaining 1's
            new_tcl(datasample(ind_zer,remo,'Replace',false)) = 1; % assign remaining ones
        else 
            new_tcl = tcl;
            ind_one = find(new_tcl == 1);
            remo = round((sum(tcl)) - sum(ziE(xi,:))); % number of overesimated 1's
            new_tcl(datasample(ind_one,remo,'Replace',false)) = 0; % assign zero randomly to equalize number of 1's
        end
     ziE2(xi,:) = new_tcl;
    end

    % points variability from the theoretical Pk distribution
    if rand < options.prob_dev
        r = binornd(1,rand,1,length(ziE2(xi,:))); % random probability 0-1
        tmpz = ziE2(xi,:) - r;
        tmpz(tmpz < 0) = 0;
        ziE2(xi,:) = tmpz;
    end

end

% now simulate sequencing noise error
if strcmp(options.noise_model,"ZINB") 
    E = random('Negative Binomial',1/options.SNR,0.5,size(P));
    X = X + E;
    %E = var(X(:))/(SNR*var(nbE(:)));
end

if options.raw_counts
    Y = X.*ziE2; % inflate with dropouts
    Y(Y<0) = 0; % change netagtive counts to 0, as were drawn from gaussian that exceeds 0 
else
    Y = X.*ziE2; % inflate with dropouts
    Y(Y<0) = 0; % change netagtive counts to 0, as were drawn from gaussian that exceeds 0
    Y = Y./sC; % going back to fold change
    Y = log10(Y); % going back to log 10 fold change
end

E = E.*ziE2;

end