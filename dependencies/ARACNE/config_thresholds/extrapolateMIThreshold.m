function [threshold coef] = extrapolateMIThreshold(data, h, n,method, varargin)

% EXTRAPOLATEMITHRESHOLD extrapolates mutual information (MI) thresholds for 
% different statistical significance (p-values) under the the null 
% hypotheis of mutually independent variables.
%
%   [THRESHOLD COEF] = EXTRAPOLATEMITHRESHOLD(DATA, H, N) computes N 
%   pairwise MI between randomly selected variables whose sample labels are
%   permutated at random. It is then used as a null distribution for 
%   assessing the statistical significance of MI computed from the same
%   dataset. By large deviation theory, the right tail of this null
%   distriubtion should decay exponentially, which can be used to
%   extrapolate to arbitrarily small p-values.
%       DATA        - a N-by-M matrix where variables are in rows and 
%                     samples in columns.
%
%       N           - number randomly selected varialbe pairs to compute 
%                     under null
%
%       H           - kernel width used in MI calculation.
% 
%       THRESHOLD   - MI threshold determined for p-values of 
%                     {5E-1, 2E-1, 1E-1, 5E-2, 2E-2 1E-2, ..., 1E-20}
%
%       COEF        - coefficient of the extrapolation function fitted to
%                     the right tail of the null distribution
%
%   [THRESHOLD COEF] = EXTRAPOLATEMITHRESHOLD(DATA,H)N,'PLOT',T), where 'T=1' 
%   displays a plot of the extrapolation; 'T=0' otherwise.  
%

graphplot = 0;

if nargin < 4
    error('Incorrect number of arguments to PVALUETABLE.');
end
if rem(nargin-4,2)~=0
    error('Incorrect number of optional arguments to PVALUETABLE.');
end

okargs = {'plot'};
for it=1:2:length(varargin)
    paraname = varargin{it};
    paraval = varargin{it+1};
    k = strmatch(lower(paraname), okargs);
    if isempty(k)
        error('Unknown parameter name: %s.', paraname);
    else
        if ~isnumeric(paraval)
            error('Display can be either 1 or off 0.');
        else
            graphplot = paraval;
        end
    end
end

% permute the genes to generate random pairs
% reset the random number generator back to the same initial state
%rand('state', 0); 
[ng,ma] = size(data);
rdata   = zeros(ng,ma);
for i=1:ng
    rdata(i,:) = data(i,randperm(ma));
end

% compute null MI 
it = 1;
null = zeros(n,1);
clear global;
while (it <= n) 
    i = round(rand*(ng-1)+1);
    j = round(rand*(ng-1)+1);
    if (i == j)
        while (j == i) 
            j=round(rand*(ng-1)+1);
        end
    end
    if strcmp(method,'fixed_bandwidth')
        null(it) =  mutualinfo_fb(rdata(i,:),rdata(j,:),h,'cleartable',0);
    end
    if strcmp(method,'variable_bandwidth')
        sigmaX = kernelDensityBandwidth(rdata(i,:));
        sigmaY = kernelDensityBandwidth(rdata(j,:));
        null(it) = mutualinfo_vb(rdata(i,:),rdata(j,:),sigmaX,sigmaY); 
    end
    if strcmp(method,'adaptive_partitioning')
        null(it) = mutualinfo_ap([rdata(i,:)' rdata(j,:)']); 
    end
    it = it+1;
end    

[c,I] = sort(null);
idx   = 100:100:n;
cumf  = zeros(length(idx),2);
for i=1:length(idx)
    cumf(i,1) = c(idx(i));
    cumf(i,2) = (n-idx(i))/n;    
end

nfit = 100;
% ids  = find(cumf(:,1)>0.03 & cumf(:,2)~=0);
foo  = find(cumf(:,2)~=0);
if length(foo)<nfit
    error('not enough points to fit.');
end
ids  = foo((end-nfit+1):end);
xx   = [ones(length(ids),1) cumf(ids,1)];
yy   = log(cumf(ids,2));
coef = regress(yy,xx);

% plot graph
if (graphplot==1)
    figure;
    subplot(2,1,1);
    [b,a] = hist(null, 100);
    b = b./n;
    plot(a, b,'LineWidth',2);
    xlabel('\it{I_0}','FontSize',20);
    ylabel('Density','FontSize',20);
    set(gca,'FontSize',20);
    
    subplot(2,1,2);
    semilogy(cumf(:,1),cumf(:,2),'b.','MarkerSize',15);
    hold on;
    yyy = 0:-1:-round(log(n));
    xxx = (yyy-coef(1))./coef(2);
    semilogy(xxx,power(exp(1),yyy),'r-','LineWidth',2);
    hold off;
    xlabel('\it{I_0}','FontSize',20);
    ylabel('\it{p}-value','FontSize',20);
    set(gca,'FontSize',20);
    pause(0.1);
    awtinvoke(get(gcf,'javaframe'),'maximize');
    print('-painters','-dmeta','-r300',['figure_mi_null_extrapolation_' num2str(size(data,2)) '.wmf']);
    close;
end

%compute p-value vs thresholds
p           = 20; 
decimal     = [5 2 1];
threshold   = zeros(p*length(decimal),1);
% null_sorted = sort(null);
% for j=1:p
%     for k = 1:length(decimal)
%         pv = decimal(k)*10^(-j);
%         % since we compute only 1e+5 mutual informations, the p-values can only
%         % be accurately estimate up to 1e-4; the rest need to extrapolate
%         if (j<5)    
%             threshold(3*(j-1)+k) = null_sorted(ceil((1-pv)*n));
%         else
%             threshold(3*(j-1)+k) = (log(pv)-coef(2))/coef(1);
%         end
%     end
% end
