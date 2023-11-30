function mi = mutualinfo_vb(xdata,ydata,sigmaX,sigmaY,xEvalPoints,yEvalPoints)
% Copyright 2008 Mukesh Bansal

%% Converting data to column format
if size(xdata,1)==1
    xdata = xdata';
end
if size(ydata,1)==1
    ydata = ydata';
end

%% Parse input parameters
if nargin < 4
    error('Incorrect number of required arguments to MUTUALINFO.');
end
if (length(xdata) ~= length(ydata))
    error('Vectors must be of the same length!');
end
if nargin ==4
    xEvalPoints = xdata;
    yEvalPoints = ydata;
end


%% Initializing parmaters
mu = [0 0];
ss = 0;
ss1 = 0;
n = length(xdata);
covarianceMatrix = [sigmaX^2 0;0 sigmaY^2];

%% Computing mutual information
dx = abs(xdata*ones(1,n) - ones(n,1)*xEvalPoints');
dy = abs(ydata*ones(1,n) - ones(n,1)*yEvalPoints');
fxy = sum(multinormal([dx dy], mu, covarianceMatrix,n)) / length(xdata);
fx = sum(normal(dx, mu(1), sigmaX)) / length(xdata);
fy = sum(normal(dy, mu(2), sigmaY)) / length(ydata);
ss = sum(log(fxy./(fx.*fy)));
minfo = max(ss/n, 0);
mi = minfo;

%% Compute normal pdf
function y = normal(dx,mu,sigma)
y = exp(-0.5 * ((dx - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);

%% Compute multivariate pdf
function y = multinormal(data,mu,covariancematrix,n)
y = exp(-.5*(((data(:,1:n) - repmat(mu(1),size(data(:,1:n))))/sqrt(covariancematrix(1,1))).^2 ...
    + ((data(:,n+1:2*n) - repmat(mu(2),size(data(:,n+1:2*n))))/sqrt(covariancematrix(2,2))).^2));
y = y/(2*pi*sqrt(covariancematrix(1,1))*sqrt(covariancematrix(2,2)));


