function export2gnuplot(file,varargin)
% function that exports vectors and variables to gnuplot tsv data file
% format.
% 
% export2gnuplot(file,variables,data)
% 
% Input arguments:
% ===============
% file: path and filenname to store data in. Path must exist. If file
%       exist the file will be appended with the new data
%
% variable input arguments:
% ========================
% variables: variable names if any. Cell array with strings
% data: 2d matrix, where each row represent a data point.
% 
% File structure:
% ==============
% 
% # Header
% # variable{1} variable{2} ... variable{n-1} variable{n}
%   data(1,1)   data(1,2)   ... data(1,n-1)   data(1,n)
%       .           .       ...     .             .
%       .           .       ...     .             .
%   data(n,1)   data(n,2)   ... data(n,n-1)   data(n,n)
% 
% 
% Header is defined by having an empty data source for the first variable data pair. 
%     

    