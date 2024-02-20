function [E, stdE] = noise(X, P, SNR, options)

    arguments
        X double % noise-free data
        P double % P matrix
        SNR double = 0.1 % value of signal-to-noise ratio
        options.method (1,1) {mustBeText} = 'SNR_L'
    end 

    if options.method == "SNR_L"
        s = svd(X);
        stdE = min(s)/(SNR*sqrt(chi2inv(1-analyse.Data.alpha,numel(P))));
        E = stdE*randn(size(P));
    elseif options.method == "SNR_var"
        stdE = sqrt(var(X(:))/SNR);
        E = stdE*randn(size(P));
    elseif options.method == "SNR_wiki"
        E = zeros(size(P));
        for i = 1:(size(X,1))
            tmpX = X(:, P(i,:)~=0);
            m = mean(tmpX(:));
            stdE = abs(m/SNR);
            E(:, P(i,:)~=0) = stdE*randn(size(tmpX));
        end
    elseif options.method == "SNR_wiki2"
        E = zeros(size(P));
        for i = 1:(size(X,1))
            tmpX = X(:, P(i,:)~=0);
            m = mean(tmpX(:));
            stdE = sqrt(abs((m^2)/SNR));
            E(:, P(i,:)~=0) = stdE*randn(size(tmpX));
        end
    elseif options.method == "SNR_cov"
        E = zeros(size(P));
        for i = 1:(size(X,1))
            tmpX = X(:, P(i,:)~=0);
            m = mean(mean(cov(tmpX)));
            stdE = sqrt(abs((m^2)/SNR));
            E(:, P(i,:)~=0) = stdE*randn(size(tmpX));
        end
    elseif options.method == "SNR_Lrepcor" % experimental code, works only for 2 replicates data
        E = zeros(size(P));
        for i = 1:(size(X,1))
            tmpX = X(:, P(i,:)~=0);
            N = size(tmpX,1);
            s = svd(tmpX);
            stdE = s(N)/(SNR*sqrt(chi2inv(1-analyse.Data.alpha,numel(P))));
            Er = stdE*randn(size(tmpX));
            Er0 = (Er(:,2)*SNR) + (sqrt(1-SNR^2))*Er(:,1); % r is replaced by SNR, so here SNR is expected correlation value 
            % Er0 will be a new variable correlated with Er(:,2) with p = SNR
            E(:, P(i,:)~=0) = [Er0, Er(:,2)];
        end
    elseif options.method == "SNR_repcorvar" % same as above
        E = zeros(size(P));
        for i = 1:(size(X,1))
            tmpX = X(:, P(i,:)~=0);
            stdE = sqrt(var(X(:))/SNR);
            Er = stdE*randn(size(tmpX));
            Er0 = (Er(:,2)*SNR) + (sqrt(1-SNR^2))*Er(:,1);
            E(:, P(i,:)~=0) = [Er0, Er(:,2)];
        end
    elseif options.method == "SNR_manual"
        stdE = SNR;
        E = stdE*randn(size(P));
    end



end




    
