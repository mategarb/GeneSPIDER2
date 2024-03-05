function varargout = BC3NET(varargin)

rawZeta = 0;
regpath='input';
net = [];

for i=1:nargin
    if isa(varargin{i},'datastruct.Dataset')
        data = varargin{i};
    elseif isa(varargin{i},'datastruct.Network')
        net = varargin{i};
    elseif isa(varargin{i},'char')
        regpath = varargin{i};
    elseif isa(varargin{i},'logical')
        rawZeta = varargin{i};
    else
%<<<<<<< HEAD
        %if length(varargin{i}) > 1
%=======
        if length(varargin{i}) > 1 | varargin{i} > 1
%>>>>>>> 2020-12-16_Deniz_genespider
            zetavec = varargin{i};
        else
            alpha = varargin{i};
        end
    end
end

if ~exist('data')
    error('%s needs a data set', mfilename)
end

% make some tmp files and directories that will allow us to feed
%<<<<<<< HEAD
% the data between matlab and R 
% this is done dynamically so that tmp files are always in pwd 
filewp = mfilename('fullpath');
filewp = filewp(1:end-6);

cwd = pwd; 

tmpdir = strcat(cwd, '/BC3NET'); 
%=======
% the data between matlab and R
% this is done dynamically so that tmp files are always in pwd
filewp = mfilename('fullpath');
filewp = filewp(1:end-6);

cwd = pwd;

tmpdir = strcat(cwd, '~/Benchmark/BC3NET');
%>>>>>>> 2020-12-16_Deniz_genespider
if ~exist(tmpdir,'dir')
    mkdir(tmpdir);
end

%<<<<<<< HEAD
tmp_file = strcat(tmpdir,'/data.csv');

csvwrite(tmp_file, data.Y);

call_on = strcat(filewp,'/bc3net.R');

% once all files are generated send the info to R to run bc3net 
system(['Rscript ' call_on, ' ', tmp_file,' ', tmpdir])

% and then read the output back in to matlab and remove the tmpfiles
BNET = readtable(strcat(tmpdir,'/bnet.csv'));
%=======
tmp_file = strcat(tmpdir,'~/Benchmark/BC3NET/data.csv');

csvwrite(tmp_file, data.Y);

call_on = strcat(filewp,'~/Benchmark/bc3net.R');

% once all files are generated send the info to R to run bc3net
system(['Rscript ' call_on, ' ', tmp_file,' ', tmpdir])

% and then read the output back in to matlab and remove the tmpfiles
BNET = readtable(strcat(tmpdir,'~/Benchmark/BC3NET/bnet.csv'));
%>>>>>>> 2020-12-16_Deniz_genespider
BNET = table2array(BNET);

system(['rm -r ',tmpdir]);

varargout{1} = BNET;

%<<<<<<< HEAD
% once this is done do a LSCO like sparsity adjustment on the data 
%=======
% once this is done do a LSCO like sparsity adjustment on the data
%>>>>>>> 2020-12-16_Deniz_genespider
if ~exist('zetavec','var') & strcmpi(regpath, 'input')
    estA = BNET;
else
    Als = BNET;
end

%<<<<<<< HEAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine how to handle zeta %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(regpath,"AZ")
%=======
if strcmpi(regpath,'AZ')
%>>>>>>> 2020-12-16_Deniz_genespider
    % test if zetavec is a array, if yes get the length
    if length(zetavec) > 1
        zetavec = length(zetavec);
    end
    % then draw a set of values that give a stepwise sparsity
    % between full and empty
    zetavec = gsUtilities.adaptive_sparsity(Als,zetavec);
%<<<<<<< HEAD
else
    if strcmpi(regpath,'full')
        zetavec = abs(Als(:)');
        zetavec = unique(zetavec);
    end

    if ~rawZeta & strcmpi(regpath,'input')
        zetaRange = [];
        estA = Als;
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
end
%=======
elseif strcmpi(regpath,'full')
    zetavec = sort(abs(Als(:)'));
    zetavec = unique(zetavec);
elseif strcmpi(regpath, 'CO')
    zetavec = sort(abs(Als(:)'));
    zetavec = unique(zetavec);
    z = randperm(length(zetavec), 18);
    zetavec = [min(zetavec) sort(zetavec(z)) max(zetavec)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine how to handle zeta %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~rawZeta & strcmpi(regpath,'input')
    zetaRange = [];
    estA = Als;
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

%>>>>>>> 2020-12-16_Deniz_genespider
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(zetavec)
    temp = find(abs(Als) <= zetavec(i));
    Atmp = Als;
    Atmp(temp) = 0;
    estA(:,:,i) = Atmp;
end
varargout{1} = estA;

if nargout > 1 & strcmpi(regpath, 'input')
    varargout{2} = zetavec;
 %   varargout{3} = zetaRange;
elseif nargout > 1 & strcmpi(regpath, 'full')
    zetavec = (zetavec-zetaRange(1))/delta;
    varargout{2} = zetavec;
  %  varargout{3} = zetaRange;
end

return
