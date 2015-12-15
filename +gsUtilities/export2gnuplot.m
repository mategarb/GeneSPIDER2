function export2gnuplot(file,varargin)
%   function for exporting vectors and variables to a gnuplot friendly tsv format.
%   export2gnuplot(file<,variables,data>)
%
%
%   Input arguments:
%   ===============
%   file: path and filename to store data in. If file
%         exist the file will be appended with the new data
%
%   variable input arguments:
%   ========================
%   variables: variable names if any else put in empty cell array {}. Cell array with strings
%   data: 2d matrix, where each row represent a data point.
%
%   File structure:
%   ==============
%
%   # Header
%   # variable{1} variable{2} ... variable{n-1} variable{n}
%     data(1,1)   data(1,2)   ... data(1,n-1)   data(1,n)
%         .           .       ...     .             .
%         .           .       ...     .             .
%     data(n-1,1) data(n-1,2) ... data(n-1,n-1) data(n-1,n)
%     data(n,1)   data(n,2)   ... data(n,n-1)   data(n,n)
%
%
%   # Next block
%
%   ==============
%   A file Header is defined by having an empty data source for the first variable data pair.
%

    nArgs = length(varargin);
    if ~logical(mod(nArgs,2)) && nArgs == 0
        tmp = dbstack;
        error([tmp.name,' needs variables (cell array)/data (2d matrix) pairs as input, see: help gsUtilities.export2gnuplot'])
    end


    [path,f,ext] = fileparts(file);
    if ~isdir(path)
        mkdir(path);
    end

    if isempty(ext)
        ext = '.tsv'
    end

    if ~exist(fullfile(path,[f,ext]),'file')
        fid = fopen(fullfile(path,[f,ext]),'w');
    else
        fid = fopen(fullfile(path,[f,ext]),'a');
        fprintf(fid,'\n\n');
    end

    if fid == -1
        error('File %s could not be opened for writing',fullfile(path,f,ext))
    end

    for i=1:2:nArgs
        variables = varargin{i};
        values = varargin{i+1};

        fprintf(fid,'#');
        if iscell(values)
            fprintf(fid,' %s\t',variables{1,:});
        elseif ~isempty(variables)
            fprintf(fid,'\t%s',variables{1,:});
        end
        fprintf(fid,'\n');
        if isempty(values)
            continue
        end

        if iscell(values) % assume labels are in the first column and values in the second column
            for j=1:size(values,1)
                for k=1:size(values,2)-1
                    fprintf(fid,'%s',values{j,k});
                end
                fprintf(fid,'\t%g',values{j,end}(1,:));
                fprintf(fid,'\n');
            end
        else
            for j=1:size(values,1)
                fprintf(fid,'\t%g',values(j,:));
                fprintf(fid,'\n');
            end
        end


        if i ~= nArgs-1
            fprintf(fid,'\n\n');
        end
    end

    fclose(fid);
