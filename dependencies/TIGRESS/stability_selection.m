function F=stability_selection(x,y,ntf,predictorTF,R,L,alpha,method)
	
% Stability selection inspired by Meinshausen & Buehlmann, 2009
%
% INPUTS:
% - x: nexp*ntg design matrix, can be normalized by column to unit variance
% - y: nexp*1 variable to predict
% - ntf: number of TFs
% - R: number of resampling runs (scalar or vector)
% - L: number of LARS steps to run
% - alpha: randomization parameter (0<alpha<=1)
%
% OUTPUT:
% - freq: a matrix of scores of size L*ntf representing the frequency of 
%   selection of each TF over the L steps.  
%
% Jean-Philippe Vert and Anne-Claire Haury, 2012 
	
Rmax=max(R);
F=cell(length(R),1);
[n p]=size(x);
halfsize= floor(n/2);
freq =zeros(L,ntf);

for i=1:Rmax
	% Randomly reweight each variable
	xs = x.* repmat(alpha + (1-alpha)*rand(1,p),n,1);
		
    % Ramdomly split the sample in two sets
	perm = randperm(n);
	i1= perm(1:halfsize);
	i2= perm((halfsize+1):n);
	
    if strcmp(method,'spams')
        % run the randomized lasso on each sample and check which variables are selected
        param.L=L;
        param.mode=2;
        param.lambda=0.00001;
        param.lambda2=0;

        % Run Lasso in fixed-steps setting
        [dum path] =mexLasso(y(i1),xs(i1,:),param);
        
        %Store
        freq(:,predictorTF)=freq(:,predictorTF) + abs(sign(path(:,2:L+1)))';
    
        % Run Lasso in fixed-steps setting
        [dum path] =mexLasso(y(i2),xs(i2,:),param);

        %Store
        freq(:,predictorTF)=freq(:,predictorTF) + abs(sign(path(:,2:L+1)))';
    elseif strcmp(method,'lars')
        path=lars(xs(i1,:), y(i1), 'lar',L);
        %Store
        freq(:,predictorTF)=freq(:,predictorTF) + abs(sign(path(2:L+1,:)));
        
        path=lars(xs(i2,:), y(i2), 'lar',L);
        %Store
        freq(:,predictorTF)=freq(:,predictorTF) + abs(sign(path(2:L+1,:)));
    end
    for r=1:length(R)
        if i==R(r)
            F{r}=freq/(2*R(r));
        end
    end
    
end