% Demo for TIGRESS

%% Please make sure to be located in TIGRESS/.

%% Add the necessary path
addpath(genpath('./'))

fprintf('This is a demo to run TIGRESS. When it pauses, you just have to hit a key or click to continue. \n \n \n')
pause

%% Load test data
fprintf('***************************************************************** \n')
fprintf('1. Load Test data: random and small dataset to test the algorithm \n')
fprintf('***************************************************************** \n \n \n')
load dataTest 
display(dataTest)
pause

%% Run TIGRESS
fprintf('****************************************** \n')
fprintf('2. Set a few options for the algorithm... \n')
fprintf('****************************************** \n \n')
R=500;
alpha=0.3;
L=3;
display(R)
display(alpha)
display(L)
pause

fprintf('****************************** \n')
fprintf('3. Now run stability selection \n')
fprintf('****************************** \n \n')
freq=tigress(dataTest,'R',R,'alpha',alpha,'L',L);
display(freq(:,:,1))
fprintf('*************************************** \n')
fprintf('Let us plot the frequencies for gene G1 \n')
fprintf('*************************************** \n \n')
pause
plot(freq(:,:,1)','linewidth',2)
legend('TF1 (= G1)','TF2','TF3','TF4')
xlabel('Number L of LARS steps')
ylabel('Frequency with which each TF is chosen')
title('Stability selection frequencies for gene 1')
pause

%% Score edges
fprintf('******************************************** \n')
fprintf('4. Score the edges using the ''area'' method \n')
fprintf('********************************************* \n \n')
scores=score_edges(freq,'method','area','L',3);
display(scores)
pause

%% Rank edges
fprintf('********************************************* \n')
fprintf('5. Finally let us rank the edges \n')
fprintf('(We can keep the first 15 edges for instance) \n')
fprintf('********************************************* \n \n')
edges=predict_network(scores,dataTest.tf_index,'cutoff',15);
display(edges)
pause 

%% Evaluate predictions
fprintf('********************************************* \n')
fprintf('6. We may now evaluate the predictions \n')
fprintf('********************************************* \n \n')
[TPR FPR PREC REC L AUROC AUPR] = GRNInferenceEvaluation(dataTest.gold_edges_num,edges);
display(AUPR)
display(AUROC)
pause
fprintf('.................................................\n')
fprintf('That''s it! Now you can try it on real data (e.g. load data1 or load ecoli) \n')
fprintf('It will be less fast...\n')


