function [A, MIfull] = clr(data, method, n, k)
%[A, MI] = clr(data, [method, n, k])
%
%returns context likelihood of mutual information between rows of 'data'
%data == genes x experiments data matrix, or, optionally, a mutual
%   information matrix (symmetric to transpose, autodetected)
%method is one of the following:
%   'rayleigh' - better convergence on smaller networks; empirically,
%   big-network fit is a little worse than stouffer
%   'beta' - same motivation as rayleigh; finite support on two sides
%   should (in theory) allow for better right tails
%   'plos' - the heuristic used in PLoS 07
%   'normal' - similar to plos, but uses unweighted Stouffer method to
%   combine z-scores
%   'stouffer' - same as 'normal', but uses weighted Stouffer method, where weights are
%   1/(column variance); may converge to good maps faster (in terms on net size) than the above,
%   but runs slower
%   'kde' - uses epanechnikov (finite support) Kernel Density estimators to fit probability
%   distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ALL METHODS USE STOUFFER METHOD TO COMBINE MARGINAL DISTRIBUTIONS
%ALL METHODS RETURN Z-SCORES INSTEAD OF P-VALUES
%FOR FALSE DISCOVERY RATE ESTIMATION, THE INCLUDED FDR.m KNOWS HOW TO WORK
%WITH Z-SCORES!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%n == number of bins (default == 10, use 5-15, or -1 for autodetect (not necessarily the best method!))
%k == spline degree (default == 3, use 2-3)
%
%Requirements: Matlab 7.0 with Statistics toolbox.  If the MI matrix is
%supplied with the 'data' argument, then Matlab 6.5 and perhaps earlier may
%also work.
%
%Additionally, it is best if the supplied mutual information matrix were of
%high-quality, e.g. not a naive partition-bin approach, but an algorithm
%utilizing some form of smoothing.  At a minimum, naive paritioning with
%Miller-Madow correction; preferably, a spline-smoothed (supplied), adaptive-binned or
%kernel density estimator approach.

[G, E] = size(data);
if G ~= E || length(find(abs(data - data') > 1e-10)) ~= 0
	fprintf('Non-square or asymmetric input matrix; assuming unprocessed data.\n');
	fprintf('Computing mutual information - this will take a while!\n');
	clear E;
	if 1 ~= exist('k')
		k = 3;
	end
	if 1 ~= exist('n')
		n = 10;
	elseif n == -1
		fprintf('Autodetecting bin count...\t');
		h = splineHistogram(data.rma, -1, k);
		l = [];
		for i = 1:length(h)
			l(i) = length(h{i});
		end
		n = floor(median(l));
		fprintf('%d\n', n);
		if n > 25
			fprintf('Automatic n seems too large (>25); edit the code to change threshold if this is fine\n');
			fprintf('Aborting\n');
			return;
		end
	end
	MI = mi(data, n, k);
	MIfull = MI;
	SE = diag(MI);
	MI = MI - diag(diag(MI));
	RawInput = 1;
else
  MIfull = data; %return this for compatibility, even though the
                 %user actually has this matrix (since it was
                 %specified as an input!)
  MI = data;
  %extract and remove entropy from the MI matrix
  SE = diag(MI); 
  MI = MI - diag(diag(MI));
  RawInput = 0;
end

if 1 ~= exist('method')
	method = 'plos';
end

if strcmp(method, 'plos')
	A = (MI - repmat(mean(MI), G, 1)) ./ repmat(std(MI), G, 1);
	A(A < 0) = 0;
	
	%square root merely adds unneded computation w/o changing rank, but it
	%helps to scale things to pseudo-z-scores
	A = sqrt(A.^2 + (A').^2);
elseif strcmp(method, 'rayleigh')
	A = zeros(G, 'single');
	for i = 1:G
		if mod(i, 300) == 0
			fprintf('%d\n', i);
		end
		row = MI(i, setdiff(1:G, i));		
		rp = raylfit(row);
		A(i, setdiff(1:G, i)) = norminv(raylcdf(row, rp));
	end
	%new way - use Stouffer method
	A = (A + A')/sqrt(2);
elseif strcmp(method, 'normal')
	A = zscore(MI);
	A = (A + A')/sqrt(2);
elseif strcmp(method, 'stouffer')
	if RawInput == 0
		fprintf('Weighted Stouffer method requires raw data for calculation of variances\n');
		return
	end
	V = var(MI);
	w = 1./V;
	W = repmat(w, G, 1);
	A = zscore(MI);
	A = (W.*A + (W.*A)')./sqrt(W.^2 + (W.^2)');
elseif strcmp(method, 'beta') %use the double-bounded beta distribution
	A = zeros(G, 'single');
	for i = 1:G
		if mod(i, 300) == 0
			fprintf('%d\n', i);
		end
		row = MI(i, setdiff(1:G, i));
		rbound = max(max(row), SE(i)); %given the nature of the estimator, some large MIs may be slightly larger than the entropy
		row = row ./ repmat(rbound + 1e-5, 1, G-1);
		rp = betafit(double(row));
		A(i, setdiff(1:G, i)) = norminv(betacdf(row, rp(1), rp(2)));
	end
	A = (A + A')/sqrt(2);
elseif strcmp(method, 'kde')
	A = zeros(G, 'single');
	for i = 1:G
		A(i, :) = norminv(ksdensity(MI(i, :), MI(i, :), 'function', 'cdf', 'kernel', 'epanechnikov'));
	end
	A = (A + A')/sqrt(2);
else
	A = [];
	fprintf('The supported methods are beta, rayleigh, plos (classic), normal (stouffer method), stouffer (weighted), and kde\n');
	return
end
%some fitting functions produce 0 for p-values, or Inf for z-scores
%set to a high number to avoid problems with Inf comparisons, etc.
A(A == Inf) = 50;
A(A == -Inf) = -50;
A = A - diag(diag(A));
