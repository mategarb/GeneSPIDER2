function data=scale_data(X)

% Center-mean and scales to unit variance by column

muo = mean(X);
stdo = std(X);
stdo2 = 1 ./ stdo;
data = (X - ones(size(X,1),1)*muo ) * diag(stdo2) ; 
