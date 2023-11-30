function freq=tigress(data,varargin)

% Runs tigress and returns the frequency matrix
%
% Syntax 1: freq=tigress(dataset) % all option values set to default
% Syntax 2: freq=tigress(dataset,'option1',option1_value,...
%              'option2',option2_value...)
%
% REQUIRED INPUT: 
% - dataset: a structure containing expdata
% 
% OPTIONAL INPUTS:
% 
% - dataset can optionally contain the following fields : 'genenames'
%   (defaults to numbers), the names of the TFs 'tf_list' (defaults to 
%   genenames) or is set with the help of tf_index, indices to find the TFs 
%   among the genes'tf_index' (defaults to all genes or is set using tflist) 
% - R: number of resamplings that should be used to run stability selection 
%   (default=1000). Note that R can be a vector of increasing values. In
%   this case, TIGRESS will return frequencies for each of these values.
% - alpha: randomization level. alpha is a scalar st 0<alpha<=1. If
%   alpha=1, no randomization is used (default=.2)
% - L: number of LARS steps that should be considered (default=5)
% - verbose: do you want the algorithm to show you the progress? Logical.
%   (default=true)
% - parallel: should MATLAB run in parallel? (default=false). Warning:
% setting this option to true requires the parallel computing toolbox.
% - LarsAlgo: which algorithm should be used for LARS? Options are: 'spams'
% (see below), 'lars' (default) or 'glmnet'.
% 
% OUTPUT: 
% - freq: a frequency matrix of size (ntf,L,ntg,[length(R)])
%
% Example:
% load data1
% freq=tigress(data1,'R',[200 500 1000],'alpha',0.1,'L',3,'verbose',false)
% 
% Note that setting otpion 'LarsAlgo' to 'spams' requires the SPAMS 
% toolbox, as implemented by Julien Mairal. 
%
% See also : tigress_full, stability_selection, score_edges, predict_network
% 
% Anne-Claire Haury, 2011

%% Parse arguments

[expdata genenames tflist tfindices L R alpha verbose LarsAlgo ...
    parallel] = checkTIGRESSargs(data,varargin,'tigress') ;

expdata = data.expdata;
[nexp ngenes]=size(expdata);
ntf=length(tfindices);

if verbose
   fprintf('Found %d genes, %d experiments, %d transcription factors\n\n',ngenes,nexp,ntf)
  fprintf('L = %d \n',L)
    fprintf('alpha = %f \n',alpha)
    fprintf('R = %d \n',2*R)
    fprintf('Algorithm : %s \n\n',LarsAlgo)
end

% Normalize each gene to 0 mean and unit variance 
expdata =scale_data(expdata);

%% Prepare frequency matrix
freq=zeros(ntf,L,ngenes,length(R));
%% Treat target genes one by one
if parallel
    F=cell(ngenes,1);
    for itarget=1:ngenes
        F{itarget}=zeros(ntf,L,length(R));	
    end
    matlabpool open;
    parfor itarget=1:ngenes
        if verbose
            fprintf('Now working on target %d out of %d \n',itarget,ngenes)
        end
        targetname = genenames(itarget);

        %Find the TF to be used for prediction (all TF except the target if the target is a TF)
        predTF = tfindices(~ismember(tflist,targetname));

        % Extract the variables for regression
        x = expdata(:,predTF);
        y = expdata(:,itarget);

        % Run stability selection
        tmp=stability_selection(x,y,ntf,predTF,R,L,alpha,LarsAlgo);
        for r=1:length(R)
            F{itarget}(:,:,r) = tmp{r}';
        end
    end
    matlabpool close
    for itarget=1:ngenes
        for r=1:length(R)
            freq(:,:,itarget,r)=F{itarget}(:,:,r);
        end
    end
    
else 
    for itarget=1:ngenes
        if verbose
            fprintf('Now working on target %d out of %d \n',itarget,ngenes)
        end
        tic
        targetname = genenames(itarget);

        %Find the TF to be used for prediction (all TF except the target if the target is a TF)
        predTF = tfindices(~ismember(tflist,targetname));

        % Extract the variables for regression
        x = expdata(:,predTF);
        y = expdata(:,itarget);

        % Run stability selection
        tmp=stability_selection(x,y,ntf,find(predTF),R,L,alpha,LarsAlgo);
        for r=1:length(R)
            freq(:,:,itarget,r) = tmp{r}';
        end
        toc
    end
end