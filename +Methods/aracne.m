function varargout = aracne(varargin)
% function estA = aracne(data <,net,alpha,regpath>)
%
%   Input Arguments: aracne(data <,net,alpha,zetavec,regpath>)
%   ================
%   data:    datastruct.Dataset
%   net:     datastruct.Network <optional>
%   alpha:   confidense level of inference, can only be a single digit \in [0,1[. default 0.01
%   zetavec: should the confidence values be truncated at some specified values.
%            Will be ignored if regpath = 'full'
%   regpath  {'input','full'}  string to determine if we should try to create a zetavec
%            for the complete regularization path dependant on method, default 'input'.
%            Where 'input' is a zetavec determined by the user
%            and 'full' gives the best estimate of the complete regularization path
%  aracnedir if the PATH to aracne is not specified beforehand a PATH to the
%            binarary can be supplied here, must be a valid path
%
%   Output Arguments: confA, zetavec
%   =================
%   confA: the confident networks as a 3d array.
%   zetavec: the complete regularization path
%   zetaRange: will return the raw zeta range scale factors
%
% For this function to work, aracne homed direcotry has to be set. This is where executables and config files can be found.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse input arguments %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawZeta = 0;
regpath = 'input';
alpha = 1;
net = [];
aracnedir = '';
for i=1:nargin
    if isa(varargin{i},'datastruct.Dataset')
        data = varargin{i};
    elseif isa(varargin{i},'datastruct.Network')
        net = varargin{i};
    elseif isa(varargin{i},'logical')
        rawZeta = varargin{i};
    elseif isa(varargin{i},'char')
        if exist(varargin{i}) == 7
            aracnedir = varargin{i};
        else
            regpath = varargin{i};
        end
    else
        if length(varargin{i}) > 1
            zetavec = varargin{i};
        else
            alpha = varargin{i};
        end
    end
end

if ~exist('data')
    error('needs a data set')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(aracnedir)
    path1 = getenv('PATH');
    tmp = strsplit(path1,':')
    if ~ismember(aracnedir,path1)
        path1 = [aracnedir ':' path1];
        setenv('PATH', path1);
    end
end
nY = response(data);

save 'in.txt' -ascii -tabs nY;
if system('aracnify.py in.txt > in.tab') ~= 0,
    error('aracnify.py could not be used')
    return
end;

% find_aracnedir = ['ARACNEDIR=$(dirname `which aracne2`)'];
% system('export LD_LIBRARY_PATH=/usr/lib/')
% system('which aracne2 > bla.txt')
% system('aracne2')
% system(find_aracnedir)
return
% cmd = ['export LD_LIBRARY_PATH=/usr/lib/;cd ',ARACNEdir,[';./aracne2 ' ...
%                     '-i $OLDPWD/in.tab -e '] num2str(alpha) ' -o ' ...
%        '$OLDPWD/out.adj;cd -'];
cmd = ['export LD_LIBRARY_PATH=/usr/lib/; aracne2 -H dirname `which aracne2`  -i in.tab -e ' num2str(alpha) ' -o out.adj;'];

if system(cmd) ~= 0,
    return,
end;
if system('dearacnify.py out.adj > out.mat.txt') ~= 0, return, end;
return

% cleanup
% if system('rm out.adj in.txt in.tab') ~= 0, return, end;
Aest = load('out.mat.txt');


if strcmpi(regpath,'full')
    zetavec = abs(Aest(:)');
    zetavec = unique(zetavec);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine how to handle zeta %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~rawZeta & strcmpi(regpath,'input')
    zetaRange = [];
    zetaRange(1) = min(abs(estA(estA~=0)))-eps;
    zetaRange(2) = max(abs(estA(estA~=0)))+10*eps;

    % Convert to interval.
    delta = zetaRange(2)-zetaRange(1);
    zetavec = zetavec*delta + zetaRange(1);

elseif ~rawZeta & strcmpi(regpath,'full')
    zetaRange(2) = max(zetavec);
    zetaRange(1) = min(zetavec);
    delta = zetaRange(2)-zetaRange(1);
elseif rawZeta & strcmpi(regpath,'full')
    zetaRange(2) = 1;
    zetaRange(1) = 0;
    delta = 1;
end

for i=1:length(zetavec)
    temp = find(abs(Aest) <= zetavec(i));
    Atmp = Aest;
    Atmp(temp) = 0;

    estA(:,:,i) = Atmp;
end

varargout{1} = estA;

if nargout > 1 & strcmpi(regpath,'input')
    varargout{2} = zetavec;
    varargout{3} = zetaRange;
elseif nargout > 1 & strcmpi(regpath,'full')
    zetavec = (zetavec-zetaRange(1))/delta;
    varargout{2} = zetavec;
    varargout{3} = zetaRange;
end

return
