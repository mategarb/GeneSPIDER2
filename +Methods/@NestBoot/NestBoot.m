classdef NestBoot

    % nestboot output properties listed here
    properties
      binary_networks
      signed_networks
      summed_networks
      minab_networks
      bin_cutoff
      accumulated_frequency
      binned_frequency
      FP_rate_cross
      support_at_FP_cross
      area_shuffled
      area_measured
    end

    methods
        function nbout = NestBoot(data,method,nest,boot,zetavec,FDR,datadir,par,cpus)
    
        % Description:
        %	data must be GeneSPIDER (GS) compiled dataset with Y and P matrices, for example synthetic sets included in GS
        %	method can be any such "Bo" method contained within GS, i.e. Bolasso, BoRNI, Bolsco, Botlsco, BoCLR, BoGenie3
        %	boot and init are inner iterations and outer runs, respecively. Best with >1000 boot and >100 init
        %	zetavec defines sparsity for many GS methods
        %	see article for FDR calculation: https://academic.oup.com/bioinformatics/article/35/6/1026/5086392
        %
        % Input arguments:
        %	data        GS dataset structure
        %	net         GS network structure
        %	method      any GS method staring with "Bo" prefix, character
        %	boot        inner bootstrap iterations (Q), integer
        %	init        outer runs (R), integer
        %	zetavec     a vector of zetas defining network sparsity
        %	FDR         False Discovery Rate threshold, integer [%]
        %	par         run NestBoot as a parallel version, true or false
        %	cpus        number of CPUs in parallel version, integer
        %
        % Output arguments:
        %	XNET        bootstrapped intersection at a given support cutoff
        %	Ssum        sum of signum function support to be sign
        %	minAb       minimum absolute value of sign support
        %	sXNET       signed network
        %	orig_index  crossing of shufled and plain data
        %	ACC         accumulated structure support
        %	FREQ        frequency of bins
        %	FP          frequency rate for crossing
        %	supp        support for crossing, i.e. (orig_index-1)/init
        %
        % Sample code:
        %	data = datastruct.Dataset.fetch('Tjarnberg-ID1005621-D20151111-N10-E30-SNR10-IDY1005621.json');
        %	net = datastruct.Network.fetch([data.network,'.json']);
        %	method = 'Bolsco';
        %	boot = 10;
        %	init = 10;
        %	zetavec = logspace(-6,0,30);
        %	FDR = 5;
        %	par = false;
        %	cpus = 2;
        %	[XNET,Ssum,minAb,sXNET,orig_index,ACC,FREQ,FP,supp] = Methods.NestBoot(data,method,boot,init,zetavec,FDR,'~/',par,cpus);
        %	m = analyse.CompareModels(net,XNET);
        %
        % Reference:
        %   Morgan, Daniel, et al. "A generalized framework for controlling FDR in
        %   gene regulatory network inference." Bioinformatics 35.6 (2019): 1026-1032.
        %
        % See also:
        %   Methods, analyse, datastruct
        %
        init = nest;
        %%%%%%%%%%ROLL-BACK CODE%%%%%%%%%%
        esta = cell(1,init);
        Afrac = cell(1,init);
        Asign_max = cell(1,init); % agnostic sign support
        Asign_frac = cell(1,init); % explicit sign support -1 for 100% negative
        Apos = cell(1,init); % positive links
        sesta = cell(1,init);
        sAfrac = cell(1,init);
        sAsign_max = cell(1,init);
        sAsign_frac = cell(1,init);
        sApos = cell(1,init);
        
        if par % run parallel version
        
             % make sure that a value is always give to CPUS
            if ~exist('cpus','var')
                cpus = 1;
            end
            % start the parallel pool if it is not running
            if isempty(gcp('nocreate'))
                parpool(cpus);
            end
            
            parfor (u = 1:init, cpus)
                [~,Asign_frac{u},~,~,~] = Methods.NestBoot.(method)(data,zetavec,boot);% infer networks for measured data  ...,boot,init); %remove "0" for bolasso as it doesnt have this function
                tmpdata = struct(data);
                for M = 1:data.M
                    shuf = randperm(data.N); % a random permutation of the integers from 1 to n without repeating elements
                    tmpdata.Y(:,M) = tmpdata.Y(shuf,M);
                    tmpdata.E(:,M) = tmpdata.E(shuf,M);
                end
                shuffledata = datastruct.Dataset(tmpdata);
                [~,sAsign_frac{u},~,~,~] = Methods.NestBoot.(method)(shuffledata,zetavec,boot);% infer networks for shuffled data  ...,boot,init);% dont need to reshuffle
            
            end
        
        % analyse_bootstrap_runs     
        XNET = zeros(data.N,data.N,length(zetavec));
        Ssum = zeros(data.N,data.N,length(zetavec));
        minAb = zeros(data.N,data.N,length(zetavec));
        sXNET = zeros(data.N,data.N,length(zetavec));
        orig_index = zeros(1,length(zetavec));
        ACC = zeros(init+1,2,length(zetavec));
        FREQ = zeros(init,2,length(zetavec));
        FP = zeros(1,length(zetavec));
        supp = zeros(1,length(zetavec));
        Area_M = zeros(init-1,length(zetavec));
        Area_S = zeros(init-1,length(zetavec));
    
        support_cross = [];
        FP_rate_cross = [];
        area_shuff_out = [];
        area_meas_out = [];
        
            parfor (z = 1:(length(zetavec)), cpus)
                [XNETa,Ssuma,minAba,sXNETa,orig_indexa,acc,FF,FPrc,suppo,ar_sh, ar_me] = NB_FDR(Asign_frac,sAsign_frac,method,data,z,init,datadir,FDR,zetavec,boot,support_cross,FP_rate_cross, area_shuff_out, area_meas_out);
                XNET(:,:,z) = XNETa;
                Ssum(:,:,z) = Ssuma;
                minAb(:,:,z) = minAba;
                sXNET(:,:,z) = sXNETa;
                orig_index(z) = orig_indexa;
                ACC(:,:,z) = acc;
                FREQ(:,:,z) = FF;
                FP(z) = FPrc;
                supp(z) = suppo;
                Area_M(:,z) = ar_me;
                Area_S(:,z) = ar_sh;
            end
        
        else %regular verson (non-parallel)
        
            for u = 1:init %outer boostrap runs - generate shufled data
            [esta{u},Afrac{u},Asign_max{u},Asign_frac{u},Apos{u}] = Methods.NestBoot.(method)(data,zetavec,boot);%,boot,init); %remove "0" for bolasso as it doesnt have this function
            tmpdata = struct(data);
                for M = 1:data.M
                shuf = randperm(data.N);
                tmpdata.Y(:,M) = tmpdata.Y(shuf,M);
                tmpdata.E(:,M) = tmpdata.E(shuf,M);
                end
            shuffledata = datastruct.Dataset(tmpdata);
            [sesta{u},sAfrac{u},sAsign_max{u},sAsign_frac{u},sApos{u}] = Methods.NestBoot.(method)(shuffledata,zetavec,boot);%,boot,init);% dont need to reshuffle
            end
        
        % analyse_bootstrap_runs     
        XNET = zeros(data.N,data.N,length(zetavec));
        Ssum = zeros(data.N,data.N,length(zetavec));
        minAb = zeros(data.N,data.N,length(zetavec));
        sXNET = zeros(data.N,data.N,length(zetavec));
        orig_index = zeros(1,length(zetavec));
        ACC = zeros(init+1,2,length(zetavec));
        FREQ = zeros(init,2,length(zetavec));
        FP = zeros(1,length(zetavec));
        supp = zeros(1,length(zetavec));
        Area_M = zeros(init-1,length(zetavec));
        Area_S = zeros(init-1,length(zetavec));
    
        support_cross = [];
        FP_rate_cross = [];
        area_shuff_out = [];
        area_meas_out = [];
    
            for z = 1:(length(zetavec))
                [XNETa,Ssuma,minAba,sXNETa,orig_indexa,acc,FF,FPrc,suppo,ar_sh, ar_me] = NB_FDR(Asign_frac,sAsign_frac,method,data,z,init,datadir,FDR,zetavec,boot,support_cross,FP_rate_cross, area_shuff_out, area_meas_out);
                XNET(:,:,z) = XNETa;
                Ssum(:,:,z) = Ssuma;
                minAb(:,:,z) = minAba;
                sXNET(:,:,z) = sXNETa;
                orig_index(z) = orig_indexa;
                ACC(:,:,z) = acc;
                FREQ(:,:,z) = FF;
                FP(z) = FPrc;
                supp(z) = suppo;
                Area_M(:,z) = ar_me;
                Area_S(:,z) = ar_sh;
            end
        end
        
        % nestboot output
        nbout.signed_networks = sXNET;
        nbout.binary_networks = XNET;
        nbout.summed_networks = Ssum;
        nbout.minab_networks = minAb;
        nbout.bin_cutoff = orig_index;
        nbout.accumulated_frequency = ACC;
        nbout.binned_frequency = FREQ;
        nbout.FP_rate_cross = FP;
        nbout.support_at_FP_cross = supp;
        nbout.area_shuffled = Area_S;
        nbout.area_measured = Area_M;
    
        end %close opening function call
    end

    methods (Static, Hidden)
        % all the nestboot methods (functions) loaded here
        varargout = lsco(varargin);
        varargout = lscon(varargin);
        varargout = tlsco(varargin);
        varargout = lasso(varargin);
        varargout = elnet(varargin);
        varargout = genie3(varargin);
        varargout = ridgeco(varargin);
        varargout = rni(varargin);
        varargout = clr(varargin);
        varargout = zscore(varargin);
    end
end % end of class 

% remaining NestBoot functions
    
    %% NB FDR
    function [XNET,Ssum,minAb,sXNET,orig_index,accumulated,binned_freq,FP_rate_cross,support_cross, area_shuff_out, area_meas_out] = NB_FDR(Afrac,shuf_Afrac,method,data,z,init,datadir,FDR,zetavec,boot,support_cross,FP_rate_cross, area_shuff_out, area_meas_out)
    
        [accumulated,~,~,~,~,freq] = accumulate(Afrac,shuf_Afrac,init,z,data); % supover,shuover,...,overlaps_shuffle
  
        bins = (0:init)/init;
        counts = zeros(length(bins)-1,size(freq,2));

        for i = 1:size(freq,2)
            counts(:,i) = [histcounts(freq(:,i),'BinEdges',bins)]';
        end
    
        binned_freq = counts./repmat(sum(counts),size(counts,1),1);
        
        %Find crossing of shuffled and plain data, based on 5% FNDR
        [orig_index, area_shuff, area_meas] = FDRcutoff(binned_freq,FDR);
        
        area_shuff_out = [area_shuff_out; area_shuff];
        area_meas_out = [area_meas_out; area_meas];

        support_cross = [support_cross; orig_index/100]; % keep it for output
    
        tmp = sum(binned_freq(orig_index:end,:), 1);
   
        FP_rate_cross = [FP_rate_cross; tmp(2)/tmp(1)]; % keep it for output
        FP_rate_cross(isnan(FP_rate_cross)) = 0;
        
        % compute and collect networks
        [XNET,Ssum,minAb,sXNET] = networks(Afrac,shuf_Afrac,boot,init,z,data,orig_index,FDR,datadir,method,zetavec);

    end
    
    %% ACCUMULATE
    function [accumulated,supover,shuover,overlaps_support,overlaps_shuffle,freq] = accumulate(booAlink,booShuffleAlink,init,z,~)
    
        estimated_support_net = [];
        estimated_shuffle_net = [];
        overlaps_support = zeros(init,size(booAlink{1},1)*size(booAlink{1},1));
        overlaps_shuffle = zeros(init,size(booShuffleAlink{1},1)*size(booShuffleAlink{1},1));

        for iter = 1:init % or 1:length(booAlink) 
            tmp1 = booAlink{iter}(:,:,z);
            
            estimated_support_net = cat(1, estimated_support_net,tmp1(:));
            overlaps_support(iter,:) = tmp1(:);
    
            tmp2 = booShuffleAlink{iter}(:,:,z);
            estimated_shuffle_net = cat(1, estimated_shuffle_net,tmp2(:));
            overlaps_shuffle(iter,:) = tmp2(:);
    
        end

        freq = [estimated_support_net,estimated_shuffle_net];
    
        supover = zeros(1,init+1);
        shuover = zeros(1,init+1);
        for k = 1:init
            tmp_intersect1 = sum(mand((overlaps_support >= k/init)));
            tmp_union1 = sum(mor((overlaps_support >= k/init)));
            supover(k+1) = tmp_intersect1/tmp_union1; 
    
            tmp_intersect2 = sum(mand((overlaps_shuffle >= k/init)));
            tmp_union2 = sum(mor((overlaps_shuffle >= k/init)));
            shuover(k+1) = tmp_intersect2/tmp_union2;
        end
    
        shuover(isnan(shuover)) = 0;
        supover(isnan(supover)) = 0;
        accumulated = [supover',shuover'];
       
    end
    
    %% FDR CUTOFF
    function [orig_index, area_shuff, area_meas] = FDRcutoff(binned_freq,FDR)
        area_shuff = zeros(1, (length(binned_freq)-1));
        area_meas = zeros(1, (length(binned_freq)-1));
        for ni = 1:(length(binned_freq)-1)
            area_shuff(ni) = trapz((length(binned_freq)-ni):length(binned_freq),binned_freq(end-ni:end,2)); % shuffled
            area_meas(ni) = trapz((length(binned_freq)-ni):length(binned_freq),binned_freq(end-ni:end,1)); % measured
        end
        % divde shuffled by measured to estimate ratio
        sm_ratio = area_shuff./area_meas; %shuffled/measured

        nlks = nnz(sm_ratio <= FDR/100); %nlinks find a cutoff below given FDR
        orig_index = length(binned_freq)-nlks;
    end
    
    %% MATRIX OR
    function out = mor(in,dim)
        if ~exist('dim','var')
            dim = 1;
        end
            ndim = length((size(in)));
        if ndim < dim
            error('Index exceeds matrix dimensions.\n\n input have no dimension %d',dim)
        end
    indim = ['in(:',repmat(',:',1,ndim-1),')'];
    indim(2+dim*2) = '1';
    tmp = eval(indim);
    indim(2+dim*2) = '2';
    out = or(tmp,eval(indim));
    indim(2+dim*2) = 'i';
        for i = 3:size(in,dim)
            out = or(out,eval(indim));
        end
    end
    
    %% MATRIX AND
    function out = mand(in,dim)
    
    in(isnan(in)) = 0;
        if ~exist('dim','var')
            dim = 1;
        end
        ndim = length((size(in)));
        if ndim < dim
            error(' Index exceeds matrix dimensions.\n\n input have no dimension %d',dim)
        end
    indim = ['in(:',repmat(',:',1,ndim-1),')'];
    indim(2+dim*2) = '1';
    tmp = eval(indim);
    indim(2+dim*2) = '2';
    out = and(tmp,eval(indim));
    indim(2+dim*2) = 'i';
        for i = 3:size(in,dim)
            out = and(out,eval(indim));
        end
    end
    
    %% NETWORKS
    function [XNET,Ssum,minAb,sXNET] = networks(Afrac,~,boot,init,z,~,orig_index,FDR,~,~,~) %data,...,datadir,method,zetavec
    
        % record of every signed direction support
        booAsign = restructurematrix(Afrac);
    
        % rules for cutoff
        %if orig_index == init
        %    cutoff = (orig_index+1)/100; % remove everything when support == init (as init defines bins)
        %else
            cutoff = orig_index/100; % otherwise use the estimated bin for which FDR is 5%
        %end
        X_network = network_at_support_level(cutoff,booAsign{z});
        
        % the bootstrapped intersection at given support cutoff
        XNET = double(mand(X_network,3));
        
        % % take sign of the sum of signum function support to be sign
        Ssum = sign(sum(booAsign{z},3)); %OK
        % take the minimum absolute value of the sign support 
        minAb = min(abs(booAsign{z}),[],3); %OK
        % sXNET = double(mand(X_network2,3));
        
        % signed network creation
        net_sum = (sum(X_network,3));
        net_sum_pos = net_sum > ((boot*init) - FDR)/boot;
        net_sum_neg = net_sum < (-((boot*init) - 1)/boot);
        net_sum_neg = net_sum_neg*-1;
        sXNET = net_sum_pos + net_sum_neg;
   
    end
    
    %% MATRIX RESTUCTURING
    function [nevomatrix] = restructurematrix(matrix) %restructures matrix by gathering outputs according to their zeta values
        nit = length(matrix); % e.g. all zetas=0 are together in one cell, normally range of zetas is in one cell
        nevomatrix = cell(nit);
        %if(size(matrix{1},3) == nit)
            
            for i = 1:size(matrix{1},3)
                for j = 1:nit
                    nevomatrix{i}(:,:,j) = matrix{j}(:,:,i);
                end
            end
    end
    
    function [freq,bins] = calc_bin_freq(matrix,init)
        [counts,bins] = histcounts(matrix(:),'BinEdges',(0:init)/init);
        freq = counts./repmat(sum(counts),size(counts,1),1);
    end
    
    %% NETWORK AT SUPPORT LEVEL
    function [cutoff_net] = network_at_support_level(support_level,support_network)
        cutoff_net = support_network;
        cutoff_net2 = support_network;
        cutoff_net(cutoff_net < support_level) = 0;
        cutoff_net2(cutoff_net2 > -(support_level)) = 0;
        cutoff_net = cutoff_net + cutoff_net2;
    end