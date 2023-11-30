function generate_mutual_threshold_configuration(data,method)
% This script generates the ARACNE configuration file - "config_threshold.txt"
%
% It requires input variable 'data', in which rows represent variables and columns samples
% and a variable 'method' that must be 'fixed_bandwidth' | 'variable_bandwidth' | 'adaptive_partitioning'
% if no method is speficied the default 'fixex_abandwidth' will be used.
% This script require the following function in the MATLAB path or in the
% current working directory: extrapolateMIThreshold.m
%
% Note: finding the MI threshold for a dataset can take some amount
% of time. In our B-cell application, we determined the kernel width for
% sample size of {100, 120, ..., 360}, which took ~6 hours to complete on
% a desktop PC with Pentimum 4 processor.
if nargin<2
   method = 'fixed_bandwidth';
end

rand('state', sum(clock));
[ng,ma] = size(data);
p       = (0.3:0.05:1)';
n       = round(ma*p);
rep     = 3;            % number of replicate for each sample size
N       = 100000;         % number of null MI to compute

alpha   = zeros(length(n),rep);
beta    = zeros(length(n),rep);
if (strcmp(method,'fixed_bandwidth')||strcmp(method,'adaptive_partitioning'))
    data = data + abs(randn(size(data))*10^-7);
end
for i=1:length(n)
    disp(['Sample Size = ' num2str(n(i))]);
    if strcmp(method,'fixed_bandwidth')
        h = getKernelWidth(n(i));
    else
        h = 0;
    end
    for j=1:rep
        disp([' Repeat ' num2str(j)]);
        idx         = randsample(ma,n(i));
        [thr coef]  = extrapolateMIThreshold(data(:,idx),h,N,method,'plot',0);
        beta(i,j)   = coef(2);
        alpha(i,j)  = coef(1);
    end
end


figure,
plot(n,mean(beta,2),'bo','LineWidth',2);
hold on;
xlabel('Sample Size','FontSize',16);
ylabel('Slope of extrapolation fit','FontSize',16);
set(gca,'FontSize',16,'LineWidth',2);

coef = polyfit(n,mean(beta,2),1);
plot(n,n*coef(1)+coef(2),'r-','LineWidth',2);
m = sprintf('%.3f', coef(1));
c = sprintf('%.2f', coef(2));
legend('Observed', ['y = ' num2str(m) 'x ' num2str(c)]);

% output the parameters to MATLAB function
out = fopen('getMIThreshold.m','w');
fprintf(out,'%s\n\n','function thr = getMIThreshold(n,p)');
fprintf(out,'%s\n',['a = ' num2str(mean(mean(alpha))) ';']);
fprintf(out,'%s\n',['b = ' num2str(coef(1)) ';']);
fprintf(out,'%s\n\n',['c = ' num2str(coef(2)) ';']);
fprintf(out,'%s\n','thr = (log(p)-a)/(b*n+c);');
fclose(out);


% output the parameters to the ARACNE configuration file
out = fopen('config_threshold.txt','w');
fprintf(out,'%g\t%g\t%g',mean(mean(alpha)),coef(2),coef(1));
fclose(out);

