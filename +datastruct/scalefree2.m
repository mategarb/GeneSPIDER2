function A = scalefree2(N, S, options)
% Create a scalefree network with N nodes and specific sparsity
% second version of scalefree function that allows creating controlable
% in/out degree closer to real power law distributions

    arguments
        N double % number of nodes in the network
        S double % sparsity of a network
        options.alpha (1,1) {mustBeNumeric} = 1.2  % it is the power value in the power law distribution
        options.pasign (1,1) {mustBeNumeric} = 0.62 % probability of acitvation sign, default based on TRRUST for h. sapiens
    end 

spar = S - 1; % due to selfloops
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
Pd = (gt/sum(gt)).^(options.alpha);
C = 1/sum(Pd); %scaling factor so it sums up to 1
Pd = C*Pd;
cp = [0; cumsum(Pd)];
ind = find(rand > cp, 1, 'last');
grich = gu(ind);

% check how many genes left that are not attached anywhere
% to ensure that all genes have at least a single link
gdif = setdiff(gns, unique(nds));
if ~isempty(gdif)
    probs2 = repelem(1, length(gdif))/length(gdif);
else
    gdif = gns;
    % draw the gene that will attach to the "rich" one
    probs2 = repelem(1, N)/N;
    probs2(gns==grich) = 0; % exclude selected "rich" gene to avoid selfloops
end

probs2 = probs2/sum(probs2); % make sure probabilites sum to 1
cp2 = [0, cumsum(probs2)];
ind2 = find(rand > cp2, 1, 'last');
gnext = gdif(ind2);

% add only if it's unique connection, i.e. not exist already
allpl = unique([net0.from+net0.to; net0.to+net0.from]);
crpl = string(grich) + string(gnext);
if(isempty(find(allpl == crpl, 1)))
    i = i + 1;
    if rand < 1 % controling in and out degree, 1 means no-change from original BA algorithm
    net0(i,:) = {grich, gnext}; % add next link to the network
    else
    net0(i,:) = {gnext, grich}; % add next link to the network
    end
end

nds0 = [net0.from; net0.to];
[gt0, ~] = groupcounts(nds0);
dgr = (sum(gt0)/N)/2;
end

% edge list to adjacency
A = zeros(N);
nr = double(erase(net0.from, 'G')');
nc = double(erase(net0.to, 'G')');

for j = 1:length(nr)

    if rand <= options.pasign % randomize the activation sign of edge based on frequency from TRRUST
        val = 1;
    else % x chance of repression
        val = -1;
    end
    A(nc(j), nr(j)) = val; % outdegree, regulators are in columns

end
A(eye(N)==1) = -1;

return