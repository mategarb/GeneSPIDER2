function options = optionParser(options,varargin)
% options = optionParser(options,varargin)
% parses input options in to a struct

if isempty(varargin)
    return
end

dbout = dbstack();
parentFunc = dbout(2).name;

optionNames = fieldnames(options);

nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
    disp(optionNames)
    error([parentFunc,' needs propertyName/propertyValue pairs with names as above')
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
    inpName = lower(pair{1}); %# make case insensitive

    if any(strmatch(inpName,optionNames))
        options.(inpName) = pair{2};
    else
        error('%s is not a recognized parameter name',inpName)
    end
end

