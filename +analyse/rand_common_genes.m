function gns = rand_common_genes(s, gc, grps, pow)
% function allowing to randomize genes with probabilities
% s - number of overalpping genes
% gc - genes degree
% grps - names of genes

gns = [];
for i = 1:s
    probs = gc/sum(gc).^(pow);
    C = 1/sum(probs); %scaling factor so it sums up to 1
    probs = C*probs;
    cprobs = [0; cumsum(probs)];
    ind = find(rand > cprobs, 1, 'last');
    gns{i} = grps(ind);
    gc = gc(find(1:length(gc)~=ind));
    grps = grps(find(1:length(grps)~=ind));
end

end

