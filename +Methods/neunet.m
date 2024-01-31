function Aest = neunet(Data)

%
% This is a wrapper function for the neural network algorihtm for
% pretrained model
%

mdl = load('+Methods/nnmdl_v1.0.mat'); % load pretrained nn model
mdl = mdl.Mdl;

Y = Data.Y;
P = Data.P;
N = size(P,1);
newY = zeros(N, N);
for i = 1:N % flatten the Y
    pns = P(i,:)~=0;
    tmpY = Y(:,pns);
    newY(:,i) = mean(tmpY,2);
end

testY = newY(:);
label = predict(mdl, testY);
Aest = reshape(label, N, N); % edge list to adjacency matrix

end
