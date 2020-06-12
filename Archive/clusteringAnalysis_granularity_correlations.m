clear all;
clc;

%% load compiled cluster data
load('cluster_results_compiled.mat'); % load matrix of compiled cluster data
clusterDataCompile = table2array(clusterDataSet); % convert back to array for functions
subjectIDlist = unique(clusterDataCompile(:,1)); % get list of unique subject IDs
N = length(subjectIDlist);

%% set parameters
threshold = .05; % minimum proportion of events represented by a given cluster
minDuration = 60; % minimum duration of seated rest (seconds)
numberSDs = 2; % number of standard deviations beyond the mean considered outlier

%% get number of clusters per subject
for i_subject = 1:length(subjectIDlist)
    clusterData = [];
    subjectID = subjectIDlist(i_subject);
    index = find(clusterDataCompile(:,1)==subjectID);
    clusterData = clusterDataCompile(index,:);
    numClusters(i_subject) = size(clusterData,1);
    numClusters_threshold(i_subject) = size(clusterData((clusterData(:,7)>threshold),:),1);
end

%% load compiled seated rest data
filename = 'sitting_lengths_0s_only.xlsx';
rawData = importdata(filename);
subjectIDs = rawData.textdata;
data = rawData.data(:,2:end);
missingData = isnan(data(1,:));
data = data(:,~missingData);
subjectIDs = subjectIDs(~missingData);
subjectIDlist2 = str2double(subjectIDs);
for i_subject = 1:length(subjectIDlist2)
    subjectData = data(:,i_subject);
    subjectData = subjectData(~isnan(subjectData));
    subjectData = subjectData(subjectData>minDuration);
    numberEvents(i_subject) = numel(subjectData);
end
for i_cell = 1:length(subjectIDlist2)
    if any(subjectIDlist == subjectIDlist2(i_cell))
        dropSubject(i_cell) = 0;
    else
        dropSubject(i_cell) = 1;
    end
end
numberEvents(dropSubject==1) = []; % remove data for subjects not used in current analysis

%% load granularity variables
load('granularity_affect_measures.mat'); % load matrix of granularity variables per subject
load('dataSet18_networkMeasures_Pearson.mat'); % load matrix of network measures per subject

%% invert zICC values to get granularity (zInv)
zInv = subjectMeasures(:,2:4)*-1;

%% put data into a table
data = horzcat(subjectIDlist, numClusters', numClusters_threshold', zInv, networkMeasures(:,4), numberEvents');
columnNames = {'PPID','numClusters',['numClusters_' num2str(threshold*100) 'pct'],'zInv_Neg','zInv_Pos','zInv_M','numCom','numEvents'};
data_Table = array2table(data,'VariableNames', columnNames);
save('clusterGranularity_data.mat','data_Table');

%% run correlations
% between overall granularity and total number of clusters
[numClust_mGran_r,numClust_mGran_p] = corr(data(:,2),data(:,6));
% between overall granularity and thresholded number of clusters
[numClustT_mGran_r,numClustT_mGran_p] = corr(data(:,3),data(:,6));

% between positive granularity and total number of clusters
[numClust_pGran_r,numClust_pGran_p] = corr(data(:,2),data(:,5));
% between positive granularity and thresholded number of clusters
[numClustT_pGran_r,numClustT_pGran_p] = corr(data(:,3),data(:,5));

% between number of communities and total number of clusters
[numClust_numCom_r,numClust_numCom_p] = corr(data(:,2),data(:,7));
% between number of communities and thresholded number of clusters
[numClustT_numCom_r,numClustT_numCom_p] = corr(data(:,3),data(:,7));

% between total number of clusters and number of seated rests
[numClust_numRests_r,numClust_numRests_p] = corr(data(:,2),data(:,8));
% between thresholded number of clusters and number of seated rests
[numClustT_numRests_r,numClustT_numRests_p] = corr(data(:,3),data(:,8));

%% create scatter plots
figure;
scatter1 = scatter(data(:,6),data(:,2),[],rgb('MediumPurple'),'filled');
set(gca,'fontsize',14)
xlim([-1.4 0]);
xlabel('emotional granularity');
ylim([0 14]);
ylabel('number of clusters');
h1 = lsline;
h1.Color = 'k';
saveas(scatter1,'numClust_mGran_scatter_plot','tiff');

figure;
scatter2 = scatter(data(:,6),data(:,3),[],rgb('SeaGreen'),'filled');
set(gca,'fontsize',14)
xlim([-1.4 0]);
xlabel('emotional granularity');
ylim([0 14]);
ylabel(['number of clusters >' num2str(threshold*100) '%']);
h1 = lsline;
h1.Color = 'k';
saveas(scatter2,'numClustT_mGran_scatter_plot','tiff');

%% check robustness to outliers
% between overall granularity and total number of clusters
m_mGran = mean(data(:,6));
sd_mGran = std(data(:,6));
threshold_mGran_low = m_mGran-(numberSDs*sd_mGran);
threshold_mGran_high = m_mGran+(numberSDs*sd_mGran);
outliers_mGran_low = data(:,6)<threshold_mGran_low;
outliers_mGran_high = data(:,6)>threshold_mGran_high;

m_numClust = mean(data(:,2));
sd_numClust = std(data(:,2));
threshold_numClust_low = m_numClust-(numberSDs*sd_numClust);
threshold_numClust_high = m_numClust+(numberSDs*sd_numClust);
outliers_numClust_low = data(:,2)<threshold_numClust_low;
outliers_numClust_high = data(:,2)>threshold_numClust_high;

outliers = outliers_mGran_low + outliers_mGran_high + outliers_numClust_low + outliers_numClust_high;
data_noOutliers = data(~outliers,:);
data_outliers = data.*outliers;
data_outliers(data_outliers(:,1)==0,:) = [];
numberOutliers = sum(outliers);
newN = size(data_noOutliers,1);

[numClust_mGran_r_NO,numClust_mGran_p_NO] = corr(data_noOutliers(:,2),data_noOutliers(:,6));

figure;
scatter3 = scatter(data_noOutliers(:,6),data_noOutliers(:,2),[],rgb('MediumBlue'),'filled');
set(gca,'fontsize',14)
xlim([-1.4 0]);
xlabel('emotional granularity');
ylim([0 14]);
ylabel('number of clusters');
h1 = lsline;
h1.Color = 'k';
hold on;
scatter(data_outliers(:,6),data_outliers(:,2),[],rgb('MediumBlue'));
hold off;
saveas(scatter3,'numClust_mGran_scatter_plot_noOutliers','tiff');

% between overall granularity and thresholded number of clusters
m_numClustT = mean(data(:,3));
sd_numClustT = std(data(:,3));
threshold_numClustT_low = m_numClustT-(numberSDs*sd_numClustT);
threshold_numClustT_high = m_numClustT+(numberSDs*sd_numClustT);
outliers_numClustT_low = data(:,3)<threshold_numClustT_low;
outliers_numClustT_high = data(:,3)>threshold_numClustT_high;

outliers_T = outliers_mGran_low + outliers_mGran_high + outliers_numClustT_low + outliers_numClustT_high;
data_noOutliers_T = data(~outliers,:);
numberOutliers_T = sum(outliers_T);
newN_T = size(data_noOutliers_T,1);

[numClustT_mGran_r_NO,numClustT_mGran_p_NO] = corr(data_noOutliers_T(:,3),data_noOutliers_T(:,6));

figure;
scatter4 = scatter(data_noOutliers_T(:,6),data_noOutliers_T(:,3),[],rgb('Magenta'),'filled');
set(gca,'fontsize',14)
xlim([-1.4 0]);
xlabel('emotional granularity');
ylim([0 14]);
ylabel(['number of clusters >' num2str(threshold*100) '%']);
h1 = lsline;
h1.Color = 'k';
saveas(scatter4,'numClustT_mGran_scatter_plot_noOutliers','tiff');

%% run multiple regressions accounting for RR, IBI, affect
% model 1: Number of clusters as a function of number of rests and overall granularity
y_numClust = data(:,2); % total number of clusters
m1_X = [ones(size(data,1),1) data(:,8) data(:,6)]; % number of rests + overall granularity
[m1_b,m1_bint,m1_r,m1_rint,m1_stats] = regress(y_numClust,m1_X);
[m1_R2change,m1_Fchange,m1_pchange,m1_df1,m1_df2] = twostep_hierarchical_regression(y_numClust,m1_X,1);

% model 2: Number of clusters (thresholded) as a function of number of rests and overall granularity
y_numClustT = data(:,3); % thresholded number of clusters
m2_X = [ones(size(data,1),1) data(:,8) data(:,6)]; % number of rests + overall granularity
[m2_b,m1_bint,m2_r,m2_rint,m2_stats] = regress(y_numClustT,m2_X);
[m2_R2change,m2_Fchange,m2_pchange,m2_df1,m2_df2] = twostep_hierarchical_regression(y_numClustT,m2_X,1);

% model 3: Number of clusters (outliers removed) as a function of number of rests and overall granularity
y_numClust = data_noOutliers(:,2); % total number of clusters
m3_X = [ones(size(data_noOutliers,1),1) data_noOutliers(:,8) data_noOutliers(:,6)]; % number of rests + overall granularity
[m3_b,m3_bint,m3_r,m3_rint,m3_stats] = regress(y_numClust,m3_X);
[m3_R2change,m3_Fchange,m3_pchange,m3_df1,m3_df2] = twostep_hierarchical_regression(y_numClust,m3_X,1);

y_numClust_S = zscore(data_noOutliers(:,2)); % z-score variables to get standardized beta as output
m3_X_S = [ones(size(data_noOutliers,1),1) zscore(data_noOutliers(:,8)) zscore(data_noOutliers(:,6))]; % number of rests + overall granularity
[m3_b_S,~,~,~,~] = regress(y_numClust_S,m3_X_S);

% model 0: Number of clusters as a function of number of rests (to get residuals for plotting)
m0_X = [ones(size(data,1),1) data(:,8)]; % number of rests
[~,~,m0_r,~,~] = regress(y_numClust,m0_X);

figure;
scatter5 = scatter(data(:,6),m0_r,[],rgb('BlueViolet'),'filled');
set(gca,'fontsize',14)
xlim([-1.4 0]);
xlabel('emotional granularity');
ylabel(['number of clusters (residuals)']);
h1 = lsline;
h1.Color = 'k';
saveas(scatter5,'m1_numClustResid_mGran_scatter_plot','tiff');