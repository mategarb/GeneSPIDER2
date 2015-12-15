function varargout = export2Cytoscape(A,genes,outf,varargin)
% export2Cytoscape(A,genes,outf{,altA})
%
% Make a Cytoscape input file from Network class.
% This will output the interaction node fields <source> <target> <strength>
% to an output file in this directory.
%
% A: interaction matrix
% genes: cell array (n x 1) with gene names
% outf: name of output text file "<outf>.tsv"; it will overwrite any
%       existing files
% altA: alternative weights for the links. (optional)
%

importnet = false;
if isa(A,'datastruct.Network')
    if ispempty(genes)
        genes = A.names;
    end
    A = A.A;
elseif isa(A,'char')
    importnet = true;
    filepath = A;
end

if ~importnet
    if length(varargin) > 0
        altA = varargin{1};
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
            if exist('altA','var')
                alt_strength = altA(tgt,src);
            end

            if strength ~= 0
                if issparse(strength)
                    strength = full(strength);
                end
                if exist('altA','var')
                    if issparse(alt_strength)
                        alt_strength = full(alt_strength);
                    end
                end


                c_src = genes{src};
                c_tgt = genes{tgt};

                if exist('altA','var')
                    fprintf(fid,'%s\t%s\t%i\t%g\t%g\t%g\n',c_src,c_tgt,sign(strength),strength,abs(strength),alt_strength);
                else
                    fprintf(fid,'%s\t%s\t%i\t%g\t%g\n',c_src,c_tgt,sign(strength),strength,abs(strength));
                end
            end
        end % target genes
    end % source genes

    % Close file
    fclose(fid);
    return
else
    fid = fopen(filepath);
    [Mdata,count] = fread(fid,Inf);
    fclose(fid);

    tab = sprintf('\t');
    lf = sprintf('\n');
    Mdata = [Mdata; lf];
    delimiter = tab;
    Mdata = char(Mdata(:)');
    matchexpr = {'\r\n' '\r' '([ \n])\1+' ' *(\n|\t) *' '^\n'};
    replexpr = {'\n' '\n' '$1' '$1' ''};
    Mdata = regexprep(Mdata,matchexpr,replexpr);
    newlines = find(Mdata == lf); % where is the newline characters
    nlines = length(newlines);

    line = Mdata(1:newlines(1)-1);
    line = strsplit(line,delimiter);

    tmpA = zeros(length(genes),length(genes),length(line(3:end)));
    ind1 = find(strcmp(genes,line{1}));
    ind2 = find(strcmp(genes,line{2}));
    for j=1:length(line(3:end))
        tmpA(ind1,ind2,j) = str2num(line{j+2});
    end
    for i=1:nlines-1
        line = Mdata(newlines(i)+1:newlines(i+1)-1);
        line = strsplit(line,delimiter);
        ind1 = find(strcmp(genes,line{1}));
        ind2 = find(strcmp(genes,line{2}));
        for j=1:length(line(3:end))
            tmpA(ind2,ind1,j) = str2num(line{j+2});
        end

    end

    varargout{1} = tmpA;
end
