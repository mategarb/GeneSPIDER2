function varargout = checkTIGRESSargs(data,varargin,caller)


% Parses input arguments to TIGRESS. This function is new to version 2.1

% Anne-Claire Haury, 2012

%% Parse

p = inputParser;   % Create an instance of the class.
p.StructExpand = true;
p.CaseSensitive = false;
switch caller
    case 'tigress'
        p.addRequired('data', @(x)isfield(x,'expdata'));
    case 'tigress_full'
        % p.addRequired('data', @isstr);
        p.addRequired('data', @(x)isfield(x,'expdata'));
        p.addParamValue('method','area',@(x)any(strcmpi(x,{'area','original'})));
        p.addParamValue('cutoff', Inf, @isfloat);
        % p.addParamValue('name_net',[data,'_TIGRESS_predictions.txt'],@isstr);
        p.addParamValue('datapath','.',@isstr);
end

p.addParamValue('L', 5, @isfloat);
p.addParamValue('R', 1000, @isfloat);
p.addParamValue('alpha', .2, @(x)x>=0 && x<=1);
p.addParamValue('verbose',true,@islogical);
p.addParamValue('LarsAlgo','lars',@(x)any(strcmpi(x,{'spams','glmnet','lars'})));
p.addParamValue('parallel',false,@islogical);
p.parse(data,varargin{:})


%% Extract arguments


L=p.Results.L;
alpha=p.Results.alpha;
verbose=p.Results.verbose;
LarsAlgo=p.Results.LarsAlgo;
parallel=p.Results.parallel;

switch caller
    case 'tigress'
        data=p.Results.data;
        expdata=data.expdata;
        

        if ~isfield(data,'genenames')
            genenames = cellstr(num2str((1:size(data.expdata,2))'));
        else
            genenames = data.genenames;
        end

        if ~isfield(data,'tflist')
            tflist = genenames;
        else
            tflist=data.tflist;
        end

        if ~isfield(data,'tf_index')
            tfindices = 1:size(data.expdata,2);
        else
            tfindices = data.tf_index;
        end

        % Make sure that tflist and tfindices coincide
        if length(tflist)<length(tfindices)
            [bla tfindices] = ismember(tflist,genenames);
        elseif length(tflist)>length(tfindices)
            tflist=genenames(tfindices);
        end
        
        R=floor(p.Results.R/2);
        
        varargout={expdata , genenames , tflist , tfindices , L , R , ...
            alpha , verbose, LarsAlgo , parallel}; 
    case 'tigress_full'

        % If tigress_full calls this function, we also need the cutoff, the method
        % and the name of the file to write the network
        cutoff = p.Results.cutoff;
        method=p.Results.method; 
        % name_net=p.Results.name_net;
        datapath = p.Results.datapath;
        R=p.Results.R;
        % varargout={data, datapath, L, R, alpha, verbose , LarsAlgo, parallel, ...
        %     method, cutoff, name_net};
        varargout={data, L, R, alpha, verbose , LarsAlgo, parallel, ...
            method, cutoff};
end




