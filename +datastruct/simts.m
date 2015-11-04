function varargout = simts(A,p,tids,stdp,doplot)
% Simulating time series response.
% [L,h] = simts(A,p,tids,stdp,doplot)
%
% Function for simulating transcript changes for a given GRN and
% perturbation vector using the basic control theory model.
% The output is a figure and a matrix L that has all the transcript
% values for all genes.
%
% A: network matrix (n x n)
% p: perturbation vector (n x 1)
% tids: number of time points
% stdp: standard deviation of input noise at each time step
%       may be a scalar or an array of size(p)
% doplot: should we plot the time series, (logical)
%
% L: output data matrix (n x tids); states at each time point
% h: figure handle


n = size(A,1);
L = zeros(n,tids);
x = zeros(n,1); % Starting state is zero
L(:,1) = x;
if isrow
    p = p';
end

tau=min(abs(eig(A)));

timestep = tau*0.1;

if ~exist('stdp')
    stdp = 0;
end

% Simulate
for t = 1:tids
    dx = A*x + p + stdp.*randn(size(p));
    x = dx*timestep + x;
    L(:,t+1) = x;
end

varargout{1} = L;

if ~exist('doplot')
    return
end

% Make legend
Genes = cell(n,1);
for g = 1:n
    Genes{g} = ['x' int2str(g)];
end

if doplot
    % Plot
    h = figure();
    timerange = (0:tids)*timestep;
    plot(timerange,L');
    xlabel('t')
    xlim([0,max(timerange)])
    legend(Genes, 'Location', 'EastOutside');
    grid on

    if nargout == 2
        varargout{1} = h;
    end
end
