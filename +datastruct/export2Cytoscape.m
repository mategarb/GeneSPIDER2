function export2Cytoscape(A,genes,outf)
% export2Cytoscape(A,genes,outf)
%
% Make a Cytoscape input file from Network class.
% This will output the interaction node fields <source> <target> <strength>
% to an output file in this directory.
%
% A: interaction matrix
% genes: cell array (n x 1) with gene names
% outf: name of output text file "<outf>.tsv"; it will overwrite any
%       existing files
%

if isa(A,'datastruct.Network')
    if ispempty(genes)
        genes = A.names;
    end
    A = A.A;
end

% Initialize
nGenes = size(A,1);
outfile = [outf, '.tsv'];
fid = fopen(outfile,'w');

if isempty(genes)
    for i=1:nGenes
        genes{i} = sprintf(['G%0',num2str(floor(log10(nGenes))+1),'d'],i);
    end
end


for src = 1:nGenes  % Loop through the source genes (columns)
    for tgt = 1:nGenes  % Loop through the target genes (rows)
        strength = A(tgt,src);
        if strength ~= 0
            if issparse(strength)
                strength = full(strength);
            end
            c_src = genes{src};
            c_tgt = genes{tgt};
            fprintf(fid,'%s\t%s\t%i\t%g\n',c_src,c_tgt,sign(strength),strength);
        end
    end % target genes
end % source genes

% Close file
fclose(fid);
