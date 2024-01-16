function A = large_scalefree(m, grn_degree, options)
% create large scal-free network to simulate big gene expression data

    arguments
        m double % m - final size of M, i.e. number of genes 
        grn_degree double % mean node degree of sub networks
        options.min_sub_grn (1,1) {mustBeNumeric} = 10 % min. number of genes for subnetwork
        options.max_sub_grn (1,1) {mustBeNumeric} = 30 % max. number of genes for subnetwork
        options.min_alpha (1,1) {mustBeNumeric} = 1.5 % min. power parameter for scale-free overlap
        options.max_alpha (1,1) {mustBeNumeric} = 2.5 % max. power parameter for scale-free overlap
        options.min_ov (1,1) {mustBeNumeric} = 1; % min. number of overlapping genes
        options.max_ov (1,1) {mustBeNumeric} = 2; % max. number of overlapping genes
    end 

% initial parameters
Ni0 = 1;
Eout = [];

while Ni0 < m

    if Ni0 == 1
    % generate initial sub GRN
    s = randi([options.min_ov,options.max_ov],1,1); % number of overlapping genes between subnetworks, here 1-3 genes overlap
    N = randi([options.min_sub_grn,options.max_sub_grn],1,1); % randomize number of nodes (genes)
    S = grn_degree; % calculate sparsity degree
    A = datastruct.scalefree2(N, S,'alpha', options.min_alpha + (options.max_alpha-options.min_alpha)*rand); % create scale-free network
    %A = datastruct.cutSym(A, 0.6, 0.398); %(A, pin, pout) in order to remove selfloops
    A = datastruct.stabilize(A,'iaa','low'); % run stabilize
    dgnl = A(eye(N)==1);
    A(eye(N)==1) = 0;
    %pN = N; % previous N, keep for the position of next sub GRN

    % now something to ensure the common genes are high regulators
    [ar, ac] = find(A ~= 0);
    ws = A(A ~= 0);
    nms = append("G",string(Ni0:(Ni0+N-1)));
    E0 = table(nms(ar)',nms(ac)',ws); % all links from the sub grn as edge list
    Eout = [Eout; E0]; % add first sub grn to the list

    % select the highest regulators
    gns = [string(E0.Var1); string(E0.Var2)];
    [gc, grps] = groupcounts(gns); % count connectivity of genes (gc) and coresponding names (grps)
    gns0 = analyse.rand_common_genes(s, gc, grps, options.min_alpha + (options.max_alpha-options.min_alpha)*rand); % select overlapping genes with probability based on their degree
    Ni0 = Ni0 + N;
    Nin = Ni0;

    elseif (Ni0 > 1) && (Ni0 < (m - options.max_sub_grn))
    % first let's calculate how much we increase our network in order to
    % break the while loop when we met the final size
    
    % generate next sub GRN that will be merged with previous
    N = randi([options.min_sub_grn,options.max_sub_grn],1,1); % randomize number of nodes (genes)
    S = grn_degree; % calculate sparsity degree
    A = datastruct.scalefree2(N, S,'alpha', options.min_alpha + (options.max_alpha-options.min_alpha)*rand); % create scale-free network
    %A = datastruct.cutSym(A, 0.6, 0.398); %(A, pin, pout) in order to remove selfloops
    A = datastruct.stabilize(A,'iaa','low'); % run stabilize
    tmpdgnl = A(eye(N)==1); % saving info about diagonal
    dgnl = [dgnl; tmpdgnl((s+1):end)];
    A(eye(N)==1) = 0;
    [ar, ac] = find(A ~= 0);
    ws = A(A ~= 0);
    nms = append("G",string(Nin:(Nin+N-1)));
    En = table(nms(ar)',nms(ac)',ws); % all links from the sub grn as edge list
    
     % select the highest regulators
    gns = [string(En.Var1); string(En.Var2)];
    [gc, grps] = groupcounts(gns); % count connectivity of genes (gc) and coresponding names (grps)
    gnsi = analyse.rand_common_genes(s, gc, grps, options.min_alpha + (options.max_alpha-options.min_alpha)*rand); % select overlapping genes with probability based on their degree
    
    % now make overlapping genes common
    for j = 1:length(gnsi)
        En.Var1 = convertStringsToChars(replace(string(En.Var1)', string(gnsi(j)), string(gns0(j)))');
        En.Var2 = convertStringsToChars(replace(string(En.Var2)', string(gnsi(j)), string(gns0(j)))');
    end

    % in the last step we select next set of overlapping genes that will be
    % used as initial set for the next loop iteration, here it must be
    % according to the whole GRN, not just the last step
    En = [Eout; En];
    gnsn = [string(En.Var1); string(En.Var2)]; % all possible genes in current network
    s = randi([options.min_ov,options.max_ov],1,1); % new number of overlapping genes
    gnsn(erase(string(gnsn), string(gnsi)) == "") = []; % remove those that were already selected
    [gc, grps] = groupcounts(gnsn); % count connectivity of genes (gc) and coresponding names (grps)
    gnsi = analyse.rand_common_genes(s, gc, grps, options.min_alpha + (options.max_alpha-options.min_alpha)*rand); % select overlapping genes with probability based on their degree
    gns0 = gnsi; % set the initial overlapping genes that will be used in next loop
    %pN = N; % previous N, keep for the position of next sub GRN
    Eout = En;
    Ni0 = length(unique([En.Var1; En.Var2]));
    Nin = Nin + N;

    elseif (Ni0 >= (m - options.max_sub_grn))
    N = m-Ni0+s; % how many genes remain to fulfuill m number of total genes
    S = grn_degree; % calculate sparsity degree
    A = datastruct.scalefree2(N, S,'alpha', options.min_alpha + (options.max_alpha-options.min_alpha)*rand); % create scale-free network
    %A = datastruct.cutSym(A, 0.6, 0.398); %(A, pin, pout) in order to remove selfloops
    A = datastruct.stabilize(A,'iaa','low'); % run stabilize
    tmpdgnl = A(eye(N)==1);
    dgnl = [dgnl; tmpdgnl((s+1):end)];
    A(eye(N)==1) = 0;
    [ar, ac] = find(A ~= 0);
    ws = A(A ~= 0);
    nms = append("G",string(Nin:(Nin+N-1)));
    En = table(nms(ar)',nms(ac)',ws); % all links from the sub grn as edge list
    
    gns = [string(En.Var1); string(En.Var2)];
    [gc, grps] = groupcounts(gns); % count connectivity of genes (gc) and coresponding names (grps)
    gnsi = analyse.rand_common_genes(s, gc, grps, options.min_alpha + (options.max_alpha-options.min_alpha)*rand); % select overlapping genes with probability based on their degree
    
    % now make overlapping genes common
    for j = 1:length(gnsi)
        En.Var1 = convertStringsToChars(replace(string(En.Var1)', string(gnsi(j)), string(gns0(j)))');
        En.Var2 = convertStringsToChars(replace(string(En.Var2)', string(gnsi(j)), string(gns0(j)))');
    end

    % in the last step we select next set of overlapping genes that will be
    % used as initial set for the next loop iteration, here it must be
    % according to the whole GRN, not just the last step

    %pN = N; % previous N, keep for the position of next sub GRN
    %Ni0 = length(unique([En.Var1; En.Var2]));
    Eout = [Eout; En];
    break
    end

end

ind1 = Eout.Var1; % go back to indices of genes instead of their names
ind2 = Eout.Var2; % go back to indices of genes instead of their names
gnso = append("G",string(1:m));
gnsu = unique([ind1; ind2]);
gnsu2 = append("G", string(sort(double(erase(gnsu, 'G')))));
for i = 1:length(gnsu2)
    ind1(ind1 == gnsu2(i)) = gnso(i);
    ind2(ind2 == gnsu2(i)) = gnso(i);
end

inds1 = double(erase(ind1, 'G'));
inds2 = double(erase(ind2, 'G'));

Mfull = zeros(m); % go back back to adjacency
for i = 1:length(inds1)
    Mfull(inds1(i), inds2(i)) = Eout.ws(i);
end
A = Mfull; % interactions from rows to columns as for scnoise it needs expression in rows in order to avoid removing them
A(eye(m)==1) = dgnl;%betarnd(5,1,1,m); % A, here 5, is the parameters of skewness towards one, i.e. higher the value, hihger the chance of having value clos to 1

end

