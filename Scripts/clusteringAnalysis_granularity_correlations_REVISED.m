% revised to incorporate outlier detection and sensitivity analyses

clear all;
clc;

%% set parameters
parameters = 1; % set to 1 to use results from person-specific DP-GMM parameters; set to 2 to use results from standard DP-GMM parameters
clusterThreshold = 0; % set to 1 to drop small clusters with weight < .05
dropOutliers = 1; % set to 1 to drop observations with Cook's D > 4/n

%% load compiled cluster data
load('cluster_results_compiled.mat'); % load matrix of compiled cluster data from analyses using person-specific DP-GMM parameters
clusterDataCompile1 = table2array(clusterDataSet); % convert back to array for functions
load('cluster_results_setNandA_compiled.mat'); % load matrix of compiled cluster data from analyses using standard DP-GMM parameters
clusterDataCompile2 = table2array(clusterDataSet_setNandA); % convert back to array for functions
subjectIDlist1 = unique(clusterDataCompile1(:,1)); % get list of unique subject IDs
subjectIDlist2 = unique(clusterDataCompile2(:,1));

% load class separation metrics
rawData = importdata('Unsupervised Clustering - subject level clustering class separation.xlsx');
bayesError = rawData.data(:,2);
klDivergence = rawData.data(:,4);

%% get number of clusters per subject
% for analyses using person-specific parameters
for i_subject = 1:length(subjectIDlist1)
    clusterData1 = [];
    subjectID1 = subjectIDlist1(i_subject);
    index = find(clusterDataCompile1(:,1)==subjectID1);
    clusterData1 = clusterDataCompile1(index,:);
    numClusters1(i_subject,1) = size(clusterData1,1);
    numClusters1_T(i_subject,1) = size(clusterData1((clusterData1(:,7)>.05),:),1);
end

% for analyses using standard parameters
for i_subject = 1:length(subjectIDlist2)
    clusterData2 = [];
    subjectID2 = subjectIDlist2(i_subject);
    index2 = find(clusterDataCompile2(:,1)==subjectID2);
    clusterData2 = clusterDataCompile2(index2,:);
    numClusters2(i_subject,1) = size(clusterData2,1);
    numClusters2_T(i_subject,1) = size(clusterData2((clusterData2(:,7)>.05),:),1);
end

%% load compiled seated rest data set
load('seatedRest_results_aggregated_wGranularity.mat'); % completeDataSet
completeData = table2array(completeDataSet); % convert back to array for functions

%% invert zICC values to get granularity (zInv)
granularity = completeData(:,19:21)*-1;

%% get number of periods of seated rest (events), mean posterior probabilities
numRests = completeData(:,2);
clusterProb = completeData(:,18);

%% put data into a table
data = horzcat(subjectIDlist1, numClusters1, granularity, numRests, clusterProb);
columnNames = {'PPID','numClusters','zInv_Neg','zInv_Pos','zInv_M','numRests','clusterProb'};
data_Table = array2table(data,'VariableNames', columnNames);
%save('clusterGranularity_data.mat','data_Table');

%% remove participants with fewer than 35 seated rest periods
index = find(numRests<35); % find participants with fewer than minimum number of rests
granularity(index,:) = []; 
numRests(index) = []; 
numClusters1(index) = []; 
numClusters1_T(index) = [];
clusterProb(index) = []; 
bayesError(index) = [];
klDivergence(index) = [];

%% get descriptive stats for class separation metrics
clusterProb_mean = mean(clusterProb);
clusterProb_sd = std(clusterProb);

bayesError_mean = mean(bayesError);
bayesError_sd = std(bayesError);

klDivergence_mean = mean(klDivergence);
klDivergence_sd = std(klDivergence);

%% set analysis
if parameters == 1 && clusterThreshold == 1
    numClusters = numClusters1_T;
elseif parameters == 2 && clusterThreshold == 0
    numClusters = numClusters2;
elseif parameters == 2 && clusterThreshold == 1
    numClusters = numClusters2_T;
else
    numClusters = numClusters1;
end

%% run zero-order correlations
% number of clusters across analyses
[a1a_a2a_r,a1a_a2a_p] = corr(numClusters1,numClusters2);
[a1b_a2b_r,a1b_a2b_p] = corr(numClusters1_T,numClusters2_T);

% number of clusters
[nClust_mGran_r,nClust_mGran_p] = corr(numClusters,granularity(:,3)); % between total number of clusters and overall granularity 
[nClust_nRests_r,nClust_nRests_p] = corr(numClusters,numRests); % between total number of clusters and number of seated rests

%% run multiple regression accounting for number of rests
% model number of clusters as a function of number of rests and overall granularity
y = numClusters; 
predictors = [numRests granularity(:,3)]; 
mdl = fitlm(predictors,y);
if dropOutliers == 1
    outliers = find(mdl.Diagnostics.CooksDistance>4/mdl.NumObservations); % use Cook's distance to identify outliers
    mdl = fitlm(predictors,y,'exclude',outliers); % re-run model excluding outliers
end
    
% z-score variables to get standardized beta as output
y_S = zscore(numClusters); 
predictors_S = [zscore(numRests) zscore(granularity(:,3))]; 
mdl_S = fitlm(predictors_S,y_S);
if dropOutliers == 1
    mdl_S = fitlm(predictors_S,y_S,'exclude',outliers);
end

% use alternative multiple regression functions to get R^2 and F stats
if dropOutliers == 1
    y_NO = y;
    y_NO(outliers) = [];
    numRests_NO = numRests;
    numRests_NO(outliers) = [];
    granularity_NO = granularity;
    granularity_NO(outliers,:) = [];
    predictors_padded_NO = [ones(length(numRests_NO),1) numRests_NO granularity_NO(:,3)]; 
    [mdl_b,mdl_bint,~,~,~] = regress(y_NO,predictors_padded_NO);
    [mdl_R2change,mdl_Fchange,mdl_pchange,mdl_df1,mdl_df2] = twostep_hierarchical_regression(y_NO,predictors_padded_NO,1);
else
    predictors_padded = [ones(length(numRests),1) numRests granularity(:,3)]; 
    [mdl_b,mdl_bint,~,~,~] = regress(y,predictors_padded);
    [mdl_R2change,mdl_Fchange,mdl_pchange,mdl_df1,mdl_df2] = twostep_hierarchical_regression(y,predictors_padded,1);
end
    
% model number of clusters as a function of number of rests (to get residuals for plotting)
predictors_R = numRests; 
mdl_R = fitlm(predictors_R,y);
y_residuals = table2array(mdl_R.Residuals(:,1));
if dropOutliers == 1
    y_residuals_NO = y_residuals;
    y_residuals_NO(outliers,:) = [];
end

%% create summary table
mdl_summary = [mdl.NumObservations; mdl_b(end); mdl_bint(end,1); mdl_bint(end,2); mdl_S.Coefficients.Estimate(end);...
    mdl.Coefficients.tStat(end); mdl.DFE; mdl_R2change; mdl_Fchange; mdl_df1; mdl_df2; mdl_pchange/2];
mdl_summary_Table = array2table(mdl_summary,'RowNames',{'N','b','b_CI_l','b_CI_u','B','t','t_df','R2','F','F_df1','F_df2','p'});

%% create scatter plot of emotional granularity against residuals
if dropOutliers == 1
    figure;
    scatter1 = scatter(granularity_NO(:,3),y_residuals_NO,[],rgb('DarkViolet'),'filled');
    set(gca,'fontsize',14)
    xlim([-1.4 0]);
    xlabel('emotional granularity');
    ylabel('number of clusters (residuals)');
    h1 = lsline;
    h1.Color = 'k';
    hold on;
    scatter2 = scatter(granularity(outliers,3),y_residuals(outliers),[],rgb('Tomato'),'filled');
    hold off;
else
    figure;
    scatter1 = scatter(granularity(:,3),y_residuals,[],rgb('DarkViolet'),'filled');
    set(gca,'fontsize',14)
    xlim([-1.4 0]);
    xlabel('emotional granularity');
    ylabel(['number of clusters (residuals)']);
    h1 = lsline;
    h1.Color = 'k';
end
    