function [E, stdE] = noise(X, P, options)

    arguments
        X double % noise-free data
        P double % P matrix
        options.SNR (1,1) {mustBeNumeric} = 0.1 % signal-to-noise ratio
        options.SNR_model (1,1) {mustBeText} = "SNR_L" % signal-to-noise ratio model, choose between SNR_L (default), SNR_var, SNR_wiki, SNR_wiki2, SNR_cov and SNR_manual
    end 

    if options.SNR_model == "SNR_L"
        s = svd(X);
        stdE = min(s)/(options.SNR*sqrt(chi2inv(1-analyse.Data.alpha,numel(P))));
        E = stdE*randn(size(P));
    elseif options.SNR_model == "SNR_vov"
        stdE = sqrt(var(X(:))/ options.SNR);
        E = stdE*randn(size(P));
    elseif options.SNR_model == "SNR_movd"
        E = zeros(size(P));
        for i = 1:(size(X,1))
            tmpX = X(:, P(i,:)~=0);
            m = mean(abs(tmpX(:)));
            stdE = abs(m/options.SNR);
            E(:, P(i,:)~=0) = stdE*randn(size(tmpX));
        end
    elseif options.SNR_model == "SNR_movd2"
        E = zeros(size(P));
        for i = 1:(size(X,1))
            tmpX = X(:, P(i,:)~=0);
            m = mean(abs(tmpX(:)));
            stdE = sqrt(abs((m^2)/options.SNR));
            E(:, P(i,:)~=0) = stdE*randn(size(tmpX));
        end
    elseif options.SNR_model == "SNR_cov"
        E = zeros(size(P));
        for i = 1:(size(X,1))
            tmpX = X(:, P(i,:)~=0);
            m = mean(mean(cov(tmpX)));
            stdE = sqrt(abs((m^2)/options.SNR));
            E(:, P(i,:)~=0) = stdE*randn(size(tmpX));
        end
    elseif options.SNR_model == "SNR_Lrepcor" % experimental code, works only for 2 replicates data
        E = zeros(size(P));
        for i = 1:(size(X,1))
            tmpX = X(:, P(i,:)~=0);
            N = size(tmpX,1);
            s = svd(tmpX);
            stdE = s(N)/(options.SNR*sqrt(chi2inv(1-analyse.Data.alpha,numel(P))));
            Er = stdE*randn(size(tmpX));
            Er0 = (Er(:,2)*options.SNR) + (sqrt(1-options.SNR^2))*Er(:,1); % r is replaced by SNR, so here SNR is expected correlation value 
            % Er0 will be a new variable correlated with Er(:,2) with p = SNR
            E(:, P(i,:)~=0) = [Er0, Er(:,2)];
        end
    elseif options.SNR_model == "SNR_repcorvar" % same as above
            E = zeros(size(P));
        for i = 1:(size(X,1))
            tmpX = X(:, P(i,:)~=0);
            stdE = sqrt(var(X(:))/SNR);
            Er = stdE*randn(size(tmpX));
            Er0 = (Er(:,2)*options.SNR) + (sqrt(1-options.SNR^2))*Er(:,1);
            E(:, P(i,:)~=0) = [Er0, Er(:,2)];
        end
    elseif options.SNR_model == "SNR_manual"
            stdE = options.SNR;
            E = stdE*randn(size(P));
    end



end




    
