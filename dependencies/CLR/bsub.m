
N = 153;
P = double([eye(N)]);
cvY = cov(aa.schurGS);
cvP = cov(P);
E = zeros(size(P));
F = zeros(size(P));
Dat(1).E = E;
Dat(1).F = F;
Dat(1).Y = aa.schurGS;
Dat(1).P = P;
reps = 1;
Yvar = (std(aa.schurGS)).^2;
gl = mean(Yvar'); % mean variance of each gene
lambda = sum((reps-1)*Yvar(:)/((reps-1)*length(Yvar(:))));
lambdas = sum(Yvar(:)')/size(Yvar(:),2);
Dat(1).lambda = lambdas;
Dat(1).cvY = cvY;
Dat(1).cvP = cvP;
data = datastruct.Dataset(Dat);
bsub_net=Methods.fcls(aa.GS,data);