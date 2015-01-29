function options = optionParser(options,varargin)
% parses input options in to a struct
% options = optionParser(options,varargin)
%
%   Input arguments
%   ===============
%   options: struct containing default options and values.
%   varargin: name value pair of input options. names must be part
%             of the options namespace.
%
%   Output arguments
%   ================
%   options: struct with translated values from input name/value pairs.
%            If input name/value pairs is empty this will return the input
%            options structure.
%


if isempty(varargin)
    return
end

dbout = dbstack();
parentFunc = dbout(2).name;

optionNames = fieldnames(options);

nArgs = length(varargin);
if round(nArgs/2) ~= nArgs/2
    disp(optionNames)
    error([parentFunc,' needs propertyName/propertyValue pairs with names as above'])
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
    inpName = lower(pair{1}); %# make case insensitive

    if any(strmatch(inpName,optionNames))
        options.(inpName) = pair{2};
    else
        error('%s is not a recognized parameter name',inpName)
    end
end
