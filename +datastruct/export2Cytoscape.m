function export2Cytoscape(A,genes,outf)
%
% Make a Cytoscape input file for Taipale's Drosophila data.
% This will output the interaction node fields <source> <target> <strength>
% to an output file in this directory.
%
% A: interaction matrix
% genes: cell array (n x 1) with gene names
% outf: base name of output text file "<outf>.txt"; it will overwrite any
%       existing files
%

if strcmp(class(A),'GeneSpider.Network')
    A = A.A;
end

% Initialize
nGenes = size(A,1);
outfile = [outf '.txt'];
fid = fopen(outfile,'w');

for src = 1:nGenes  % Loop through the source genes (columns)
    for tgt = 1:nGenes  % Loop through the target genes (rows)
        strength = A(tgt,src);
        if strength ~= 0
            if issparse(strength)
                strength = full(strength);
            end
            c_src = char(genes(src));
            c_tgt = char(genes(tgt));
            fprintf(fid,'%s\t%s\t%i\t%d\n',c_src,c_tgt,sign(strength),strength);
        end
    end % target genes
end % source genes

% Close file
fclose(fid);
