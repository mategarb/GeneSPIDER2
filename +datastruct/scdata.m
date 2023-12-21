function [Y, Ed, Eg] = scdata(A, P, options)
% Create a noise matrix E based on log normal and zinb distributions
% with preferential attachment

    % argument list with defaults
    arguments
        A double % GRN
        P double % perturbation matrix
        options.SNR (1,1) {mustBeNumeric} = 0.1 % signal-to-noise ratio
        options.raw_counts (1,1) {mustBeNumericOrLogical} = true % if true, raw counts are outputed, if false fold-change
        options.left_tail (1,1) {mustBeNumeric} = 1  % power of left tail in Pk distribution (>1 makes left tail longer)
        options.right_tail (1,1) {mustBeNumeric} = 2 % power of right tail in Pk distribution (>1 makes right tail longer)
        options.negbin_prob (1,1) {mustBeInRange(options.negbin_prob,0,1)} = 0.5 % mean for negative binomial used to simulate pseudo control
        options.disper (1,1) {mustBeNumeric} = 1 % dispersion parameter in Pk (Svensson et al.)
        options.n_clusts (1,1) {mustBeNumeric} = 5 % theoretical number of clusters
        options.prob_clusts {mustBeInRange(options.prob_clusts,0,1)} = 0.5 % probability of gene belonging to a cluster
        options.prob_dev {mustBeInRange(options.prob_dev,0,1)} = 0.3 % probability of gene being deviated from Pk=0  
        options.dev {mustBeInRange(options.dev,0,1)} = 0.2 % deviation from Pk=0, 1 means the strongest possible but still random
    end  


Net = datastruct.Network(A', 'scnet');
X = Net.G*P; %corresponding response Y, G is the static gain matrix (inverse of A (network matrix))


psc = nbinrnd(1,options.negbin_prob,[1, size(P,1)]);  % pseduocontrol mean from negative binomial
psc = psc + rand(1, length(psc)); % add some 0-1 rand that will make it behave as non-discrete mean
psc = sort(psc,'ascend');
psc = psc.^options.right_tail; % to make the tail greater and leave 0s and 1s as 0s and 1s


s = svd(X);
stdE = s(size(P,1))/(options.SNR*sqrt(chi2inv(1-analyse.Data.alpha,numel(P))));
Eg = stdE*randn(size(P)); % noise matrix
X = X + Eg;

X2 = ((10.^(X))); % reverse (pseudo)log
vX = var(X,0,2); % take variance to relate it with mean
[varx, ~] = sort(vX);
[~, idx] = ismember(vX, varx);
psc2 = psc(idx); % mean expression sorted according to gene variance

sC = zeros(size(P));
for i = 1:length(psc2)
        sC(i,:) = poissrnd(psc2(i), [1, size(P,2)]); % some pseudo control values to multiply
end

X3 = round(X2.*sC); % multiply to get pseudo counts

ziE = ones(size(P));
me = abs(mean(X3,2));
dspr = options.disper;
Pk0 = ((dspr.^-1)./(me + dspr.^-1)).^(dspr.^-1); % probability of dropouts where k = 0, from the paper
Pk0(Pk0 > 1) = 1; % because some probs are slighyly greater than 1 (rounding error most probably)

for i2 = 1:length(Pk0)
    if rand < 0.3
        if rand < 0.5
            if ~isempty(nonzeros(X3(i2,:)))
                X3(i2,:) = X3(i2,:) + (round(betarnd(options.dev,1,1,length(X3(i2,:)))*rand*max(nonzeros(X3(i2,:))))); % deviation from mean expression
            else
                X3(i2,:) = X3(i2,:) + round(betarnd(options.dev,1,1,length(X3(i2,:)))*rand); % deviation from mean expression
            end
        else
            Pk0(i2) = Pk0(i2) + betarnd(options.dev,1,1,1)*rand*options.dev; % deviation from dropout rate
            if Pk0(i2) > 1 % in case we go too far
                Pk0(i2) = 1;
            end
        end
    end
   ziE(i2,:) = random('Binomial', 1, (1-Pk0(i2))^options.left_tail, [1, size(P,2)]); % zero inflation, p is the proportion of 0s
end

% simulating UMAP clusters for random genes
clst = options.n_clusts;

tots = size(ziE,2);
int1 = round(size(ziE,2)/clst);
for it = 1:clst
    if tots > int1
        invs(it) = int1;
        tots = tots - int1;
    else
        invs(it) = tots;
    end
end

if sum(invs) < size(ziE,2)
    invs(end) = invs(end) + (size(ziE,2) - sum(invs));
end

Ab = {};
for z = 1:clst
    Ab{z} = ones(invs(z), 1);
end
bdg = blkdiag(Ab{:})';


% all possible theoretical clusters
npos = nchoosek(1:clst,2);
bdg2 = zeros(size(npos,1), size(ziE,2));
for u = 1:size(npos,1)
    bdg2(u,:) = sum(bdg(npos(u,:),:));
end

bdg3 = [bdg; bdg2];
ziE2 = ziE;
for xi = 1:size(ziE,1)
    rac = rand;
    if rac < options.prob_clusts
     tcl = bdg3(randi([1,size(bdg3,1)],1),:);
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

end

% now simulate sequencing noise error
%if strcmp(options.noise_model,"ZINB") 
%    E = random('Negative Binomial',1/options.SNR,options.SNR,size(P));
%    X = X + E;
    %E = var(X(:))/(SNR*var(nbE(:)));
%end

if options.raw_counts
    Y = X3.*ziE2; % inflate with dropouts
else
    
    Y = X3.*ziE2;
    Y = Y./sC; % going back to fold change
    Y(ismissing(Y)) = 0; % NAs from 0 division
    Y(Y~=0) = Y(Y~=0) + eps; % add small value to preserve log10 links with 1 FC (from rounding error)
    Y = log10(Y); % going back to log 10 fold change
    Y(isinf(Y)) = 0; % infs from log10(1)

end

ziE3 = ziE2;
ziE3(X3 == 0) = 0;
Ed = ziE3; % final dropouts

end