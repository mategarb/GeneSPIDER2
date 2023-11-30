function h = kernelDensityBandwidth(data,varargin)
% Copyright 2008 Mukesh Bansal
%% Set Default Parameters
type = 'std'; % Median absolute deviation
prop = 1.06;                 % Gaussian
dim = 1;       % Dimension of data
%% Parse input parameters

if nargin ==0
    error('Minimum number of arguments not found');
else
    if rem(nargin-1,2)~=0
        error('Incorrect number of arguments');
    end
    okargs = {'type'};
    for it=1:2:length(varargin)
        pname = varargin{it};
        pval = varargin{it+1};
        k = strmatch(lower(pname),okargs);
        if isempty(k)
            error('Unknown parameter name: %s.', pname);
        elseif length(k)>1
            error('Ambiguous parameter name: %s', pname);
        else
            switch(k)
                case 1  % setting type of standard deviation
                    type = pval;
            end
        end
    end
end

%% calculating standard deviation
if isequal(type,'mad')
    med = median(data);
    sig = median(abs(data-med)) / 0.6745;
end

if isequal(type,'std')
    sig = std(data);            % estimate sigma (standard)
    iqrSig = .7413*iqr(data);     % find interquartile range sigma est.
    if (max(iqrSig)==0)
        iqrSig=sig;
    end;
    sig = min(sig,iqrSig);
end
%% Computing bandwidth
n = length(data);
h = prop * sig * n^(-1/(4+dim));

%% Computing interquartile range
function x = iqr_range(y)
% compute 25th percentile (first quartile)
Q1 = median(y(find(y<median(y))));

% compute 75th percentile (third quartile)
Q3 = median(y(find(y>median(y))));

% compute Interquartile Range (IQR)
x = Q3-Q1;
