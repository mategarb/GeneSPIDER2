function demo(mode)
%demo(mode)
%
%mode:
%  'regulon'
%    Runs a comparison of all CLR methods on Regulon DB sized network (not
%    as good as using the full dataset); takes < 30 minutes
%  'plos'
%    Replicates the steps to build the full network of the type used in the
%    PLoS paper (but requires M3D build 3 dataset!!!!)  Results will differ
%    slightly from the PLoS dataset.  This is how you build a large network
%    using clr.
%  'misc'
%    Demoes various miscellaneous functions in this distribution.
%
%This demo must be run from the Code/ directory

%load the compendium
load E_coli_v3_Build_3.mat
load reg_b3.mat
load tfs_b3.mat

if strmatch('regulon', mode, 'exact')
	%compare performance of different variants of CLR on regulon-sized networks
	[Astouffer, MI] = clr(data.rma(reg_b3.ps_idx, :), 'stouffer', 10, 3);
	Anormal = clr(MI, 'normal');
	Arayleigh = clr(MI, 'rayleigh');
	Abeta = clr(MI, 'beta');
	Akde = clr(MI, 'kde');
	Aplos = clr(MI, 'plos');
	
	%remove non-TF data
	zidx = setdiff(1:length(reg_b3.ps_idx), reg_b3.tfidxTest);
	Anormal(:, zidx) = 0;
	Arayleigh(:, zidx) = 0;
	Abeta(:, zidx) = 0;
	Astouffer(:, zidx) = 0;
	Akde(:, zidx) = 0;
	Aplos(:, zidx) = 0;

	[prec, sens, spec, pval] = matrixPvalue(Aplos, Akde, Arayleigh, Abeta, Anormal, Astouffer, reg_b3.Atest, length(reg_b3.tfidxTest), 95, 0.01, 100, 10);
	figure
	plot(sens*100', prec*100', 'o');
	legend('PLoS method', 'KDE (epan)', 'Rayleigh', 'Beta', 'Stouffer gauss (unweighted)', 'Stouffer gauss (weighted)');
elseif (strmatch('plos', mode, 'exact'))
	%compute fast & coarse mutual information and a rayleigh map based on it using just the genes with b-numbers this step takes 20-30 minutes on a fast machine
	gidx = strmatch('CDS', data.probe_set_type');
	[A, MI] = clr(data.rma(gidx, :), 'normal', 7, 3);

	%zero out non-tfs
	zidx = setdiff(1:length(gidx), tfs_b3.cds_idx);
	A(:, zidx) = 0;

	%find a threshold at which the generated CLR map is 60% accurate (this works only in E. coli):
	threshold = findThreshold(A, 60);
	thresholdFDR = FDR(A, 0.05); %set a threshold based on false discovery
	%rate
	
	fprintf('60 percent precision over Regulon corresponds to a threshold of  percentf\n', threshold);
	fprintf('0.05 False Discovery Rate (FDR) corresponds to %f\n', thresholdFDR);
	fprintf('Using the Regulon threshold for this demo\n');
	
	%and the following would build the same map with the automatically determined
	%bin count
	%[A, MI] = clr(data.rma(data.gidx, :), 'rayleigh', -1, 3);
	
	%if you had a mutual information table already, you could speed things up by giving it as an argument, e.g.:
	%A = clr(MI, 'rayleigh');
	
	%make a small, 'genes-only' compendium for the purpose of drawing the
	%map (compendium is used by the mapping function to calculate correlation signs)
	c.rma = single(data.rma(gidx, :));
	c.genes = data.genes(gidx);

	mapGenes({'lexA'}, c, A, threshold, 1, 1, 'LexA map', 'lexA.ps2')
	fprintf('The map is ready; find the ''lexA.ps2'' file and open it using Adobe Illustrator or another program\n');
elseif (strmatch('misc', mode, 'exact'))
	fprintf('Illustrating a few dataset-related tasks\n');
	%quick tasks
	
	%perform enrichment on a set of genes; FDR rate of 0.05; minimum GO
	%term depth of 3
	%this will take 2-3 minutes
	[sigTerms, sigTermIds, sigGenes, sigPvalues] = enrich({'lexA', 'recA', 'dinB', 'ruvC'}, 0.05, 3)
	
	%show basR vs lexA given basS (blue == basS < median; red == basS >
	%median)
	[c, X, Y] = geneVsGene({'basR', 'lexA'}, data, 'basS');
	
	%enrich a couple of genes for transcription factors
	%The following would enrich only using Regulon:
	%[tfs, foundGenes, pvalues] = tfEnrichEColi({'recA', 'dinB'}, -1);
	%And the following will enrich using regulon and the z-map at 60%
	%estimated precision
	[tfs, foundGenes, pvalues] = tfEnrichEColi({'recA', 'dinB'});
	
	%cluster chips to see which conditions are similar (note: there are
	%a lot of chips here - this step will be slow) (another note: this
	%can work for any compendium, but it should be similarly formatted):
	buildChipTree(data);

	%slow tasks - don't run except to step through!
	fprintf('This will be done on a subsample of genes (regulon-sized set).  This is coarse but (relatively) quick compared to the full dataset.\n');
	fprintf('This will take a WEEK OR LONGER with default settings - for reference only!!!\n');
	fprintf('Press any key to continue or CTRL-C to quit\n');
	pause
	%find a good choice of 60 informationally distant chips from 20 permutations:
	[bestChipIdx, prec, sens, setSize] = findBestChips(data.rma(reg_b3.ps_idx, :), 60, 20);
	
	%compare chip selection strategies (see plot in the paper) - from 20 to
	%80 chips in step 5, repeat 10 times to estimate standard deviation:
	
	%REDUCED SETTINGS
	%[prec, sens, setSize] = netSubsample(data.rma(reg_b3.ps_idx, :), 2,
	%'rayleigh', 20, 5, 25);
	[prec, sens, setSize] = netSubsample(data.rma(reg_b3.ps_idx, :), 10, 'rayleigh', 20, 5, 80);
	
else
	fprintf('Unrecognized demo mode: %s\n', mode);
end

