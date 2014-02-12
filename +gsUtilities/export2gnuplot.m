function export2gnuplot(file,varargin)
%   function that exports vectors and variables to gnuplot tsv data file
%   format.
%
%   export2gnuplot(file,variables,data)
%
%   Input arguments:
%   ===============
%   file: path and filenname to store data in. If file
%         exist the file will be appended with the new data
%
%   variable input arguments:
%   ========================
%   variables: variable names if any. Cell array with strings
%   data: 2d matrix, where each row represent a data point.
%
%   File structure:
%   ==============
%
%   # Header
%   # variable{1} variable{2} ... variable{n-1} variable{n}
%   data(1,1)     data(1,2)   ... data(1,n-1)   data(1,n)
%       .             .       ...     .             .
%       .             .       ...     .             .
%   data(n,1)     data(n,2)   ... data(n,n-1)   data(n,n)
%
%
%   # Next block
%
%   A file Header is defined by having an empty data source for the first variable data pair.

    nArgs = length(varargin);
    if round(nArgs/2) ~= nArgs/2
        disp(optionNames)
        error([parentFunc,' needs propertyName/propertyValue pairs with names as above'])
    end

    [path,f,ext] = fileparts(file);
    if ~isdir(path)
        mkdir(path);
    end

    if ~exist(fullfile(path,f,ext),'file')
        fid = fopen(fullfile(path,f,ext),'w');
    else
        fid = fopen(fullfile(path,f,ext),'a');
        fpintf(fid,'\n\n')
    end

    if fid == -1
        error('File could not be opened for writing')
    end

    for i=1:2:nArgs
        variables = varargin{i};
        values = varargin{i=1};

        fprintf(fid,'#');
        fprintf(fid,'\t%s',variables{:});
        fprintf(fid,'\n');
        if ismepty(values)
            continue
        end
        
        for j=size(values,1)
            fprintf(fid,'\t%g',values(j,:));
            fprintf(fid,'\n');
        end
    end

    fclose(fid);