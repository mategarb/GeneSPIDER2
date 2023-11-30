% function [edges scores freq]=tigress_full(dataname,varargin)
function [edges scores freq]=tigress_full(data,varargin)
% Complete TIGRESS GRN inference method.
% Runs the TIGRESS gene regulatory network inference method given 
% expression data and a list of transcription factors. Options can be 
% included (see details below). It outputs a text file containing a 
% ranked list of edges with associated probabilities.
% For more information, please see :
% Haury et al.: 'TIGRESS: Trustful Inference of Gene REgulation using 
% Stability Selection',2012 
%
% Syntax 1: tigress_full(expression,tflist) % options to default
% Syntax 2: tigress_full(expression, tflist,'option',option_value,...)
% 
% REQUIRED INPUTS: 
%  - dataname: see read_data for info on this input
%
% OPTIONAL INPUTS:
%  - R: number of resampling runs (default 1000)
%  - L: number of LARS steps at each iteration (default:5)
%  - alpha: randomization parameter alpha (in ]0,1]) (def:0.2)
%  - method: either 'original' or 'area' (default: 'area')
%  - cutoff: number of edges to predict (default: Inf)
%  - verbose: do you want the algorithm to show you the progress? 
%            (default=true)
%  - name_net: path/name.ext of file to write network into.
%              (Default='./edges.txt')
%  - datapath: where to find the data (defaults to working directory)
%
% OUTPUT:
% A file name 'edges.txt' that writes itself in the working directory.
% It contains 3 columns:
% - 1st column: transcription factors
% - 2nd column: target genes
% - 3rd column: probability of edge existence
% Note that the file is output such that the probabilities are ranked
% decreasingly.
% 
% Example:
% tigress_full(expression_file,tflist_file,'L',3,'R',500,'verbose',false)
%
% See also: tigress.m, score_edges.m, predict_network.m
%
% Anne-Claire Haury, 2012

%% Parse arguments
% [dataname, datapath, L, R, alpha, verbose , LarsAlgo, parallel, method, ...
%     cut, name_net] = checkTIGRESSargs(dataname,varargin,'tigress_full');

% [data2, L, R, alpha, verbose , LarsAlgo, parallel, method, ...
%     cut, name_net] = checkTIGRESSargs(data,varargin,'tigress_full');

[data2, L, R, alpha, verbose , LarsAlgo, parallel, method, ...
    cut] = checkTIGRESSargs(data,varargin,'tigress_full');


% [expdata genenames tflist tfindices L R alpha verbose LarsAlgo ...
    % parallel] = checkTIGRESSargs(data,varargin,'tigress') ;


%% Read inputs
data2 = data.expdata;
% data2=read_data(datapath,dataname);

%% Get frequency matrix F
freq=tigress(data2,'L',L,'R',R,'alpha',alpha,'verbose',verbose,...
    'LarsAlgo',LarsAlgo,'parallel',parallel);

%% Get scores
scores=score_edges(freq,'method',method,'L',L);

%% Write edges
edges = predict_network(scores,data2.tf_index,'genenames',data2.genenames,...
    'cutoff',cut,'name_net',name_net);