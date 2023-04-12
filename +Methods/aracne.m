function varargout = aracne(varargin)
% function estA = aracne(data <,net,DPI,regpath>)
%
%   Input Arguments: aracne(data <,net,alpha,zetavec,regpath>)
%   ================
%   data:    datastruct.Dataset
%   net:     datastruct.Network <optional>
%   DPI:     DPI tolerance, default: 1, range \in [0,1[.
%   zetavec: should the MI values be truncated at some specified values.
%            Will be ignored if regpath = 'full'
%   regpath  {'input','full'}  string to determine if we should try to create a zetavec
%            for the complete regularization path dependant on method, default 'input'.
%            Where 'input' is a zetavec determined by the user
%            and 'full' gives the best estimate of the complete regularization path
%  aracnedir {Not implemented yet} if the PATH to aracne is not specified beforehand a PATH to the
%            binary can be supplied here, must be a valid path.
%
%   Output Arguments: confA, zetavec
%   =================
%   confA: the confident networks as a 3d array.
%   zetavec: the complete regularization path
%   zetaRange: will return the raw zeta range scale factors
%
% For this function to work, aracne home directory has to be set beforehand in the PATH variable.
% This is where executables and config files for aracne2 can be found.
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

f1='./in.tab';
t{1} = 'gene';
for i=2:data.M+1, t{i} = strcat('S',num2str(i-1)); end
delete(f1);
fix_aracne_input(f1,data);

system(['sed -i -e "s/# //" -e "s/\t$//" ', f1]);
cmd = ['export LD_LIBRARY_PATH=/usr/lib/; aracne2 -H $(dirname `which aracne2`)  -i in.tab -e ' num2str(alpha) ' -o out.adj;'];

if system(cmd) ~= 0
    return,
end

Aest = fix_aracne_output('./out.adj',data.names);

% cleanup
if system('rm out.adj in.tab') ~= 0, return, end


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

function adj = fix_aracne_output(filename,genes)

rawtxt = fileread(filename);

txtsplit = strsplit(rawtxt,'\n');

header = strncmpi('>',txtsplit,1);
txtsplit = txtsplit(~header);

adj = zeros(length(genes),length(genes));
for l=1:length(txtsplit)
    linesplit = strsplit(txtsplit{l},'\t');
    ind1 = find(strcmp(genes,linesplit{1}));
    for j=2:2:length(linesplit)
        ind2 = find(strcmp(genes,linesplit{j}));
        val = linesplit{j+1};
        adj(ind1,ind2) = str2num(val);
    end

end


function fix_aracne_input(filename,data)

nY = response(data);
outs = struct();
outs(1).gene = data.names';
for i=1:data.M
    outs(1).(['S' num2str(i)]) = nY(:,i);
end
outs = struct2dataset(outs);

export(outs,'file',filename,'Delimiter','\t')
