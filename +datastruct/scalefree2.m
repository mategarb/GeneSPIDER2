function A = scalefree2(N, n, pow, varargin)
% Create a scalefree network with N nodes and specific sparseness
% with preferential attachment
%
%
% N:    number of nodes
% n:    wished average number of links per node if n>=1, else relative sparsity of total possible links.
% pin:  The probability of adding an incoming edge to a node
% pout: pin+pout is the probabiltiy to add an outgoing edge to a given gene
%       pin and pout togheter controls the in/outgoing edges degree
%       in the powerlaw trying to ensure that the powerlaw is
%       primarily on outgoing edges. 
% seed: matrix to work as a seed, assumed to have size < N
%
% A:    undirected scalefree network matrix

% rng('shuffle');


rank_check = true;
pin = NaN;
pout = NaN;

if ~isempty(varargin)
    for i=1:length(varargin)
        if isa(varargin{i},'double')
            if isnan(pin)
                pin = varargin{i};
            elseif isnan(pout)
                pout = varargin{i};
            end
        elseif isa(varargin{i},'logical')
            rank_check = varargin{i};
        else
            seed = varargin{i};
        end
    end
end

%% parameters
%clear all
%spar = 3;
%N = 30;
%pow = 1.5;
%n=3;
spar = n;
gns = cellstr(num2str((1:N)', 'G%d'))';

lnk0 = gns(randi(N, 1, 2)); % draw two genes that will be a first seed link

varTypes = ["string","string"];
varNames = ["from","to"];
net0 = table('Size',[1 2],'VariableTypes',varTypes,'VariableNames',varNames);
net0(1,:) = lnk0; % add first link to the network

dgr = 0;
i = 1;

while dgr < spar

% draw the "rich" gene where the next will attach
nds = [net0.from; net0.to];
[gt, gu] = groupcounts(nds); %degree, gene name
Pd = (gt/sum(gt)).^(pow);
C = 1/sum(Pd); %scaling factor so it sums up to 1
Pd = C*Pd;
cp = [0; cumsum(Pd)];

ind = find(rand > cp, 1, 'last');
grich = gu(ind);


% check how many genes left that are not attached anywhere
% to ensure that all genes have at least a single link
gdif = setdiff(gns, unique(nds));
if length(gdif) ~= 0
    probs2 = repelem(1, length(gdif))/length(gdif);
else
    gdif = gns;
    % draw the gene that will attach to the "rich" one
    probs2 = repelem(1, N)/N;
    probs2(find(gns==grich)) = 0; % exclude selected "rich" gene to avoid selfloops
end

probs2 = probs2/sum(probs2); % make sure probabilites sum to 1
cp2 = [0, cumsum(probs2)];
ind2 = find(rand > cp2, 1, 'last');
gnext = gdif(ind2);

% add only if it's unique connection, i.e. not exist already
allpl = unique([net0.from+net0.to; net0.to+net0.from]);
crpl = string(grich) + string(gnext);
if(length(find(allpl == crpl)) == 0)
    i = i + 1;
    if rand < 1 % controling in and out degree, 1 means no-change from original BA algorithm
    net0(i,:) = {grich, gnext}; % add next link to the network
    else
    net0(i,:) = {gnext, grich}; % add next link to the network
    end
end

nds0 = [net0.from; net0.to];
[gt0, gr0] = groupcounts(nds0);
dgr = sum(gt0)/N;
end
% length(unique(nds0))
% edge list to adjacency
A = zeros(N);
nr = double(erase(net0.from, 'G')');
nc = double(erase(net0.to, 'G')');

for j = 1:length(nr)

    if rand <= 0.38 % randomize the sign of edge based on frequency from TRRUST
        val = -1;
    else % 0.62
        val = 1;
    end
    A(nr(j), nc(j)) = val; % outdegree

end
%A(eye(N)==1) = -betarnd(5,1,1,N); % A, here 5, is the parameters of skewness towards one, i.e. higher the value, hihger the chance of having value clos to 1
A(eye(N)==1) = -1;

%disp(length(find(real(eig(A)) > 0)))
return