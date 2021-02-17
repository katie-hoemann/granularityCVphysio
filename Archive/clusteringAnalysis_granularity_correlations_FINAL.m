clear all;
clc;

%% load compiled cluster data
load('cluster_results_compiled.mat'); % load matrix of compiled cluster data
clusterDataCompile = table2array(clusterDataSet); % convert back to array for functions
subjectIDlist = unique(clusterDataCompile(:,1)); % get list of unique subject IDs
N = length(subjectIDlist);

%% get number of clusters per subject
for i_subject = 1:length(subjectIDlist)
    clusterData = [];
    subjectID = subjectIDlist(i_subject);
    index = find(clusterDataCompile(:,1)==subjectID);
    clusterData = clusterDataCompile(index,:);
    numClusters(i_subject) = size(clusterData,1);
end

%% load compiled seated rest data set
load('seatedRest_results_aggregated_wGranularity.mat'); % completeDataSet
completeData = table2array(completeDataSet); % convert back to array for functions

%% invert zICC values to get granularity (zInv)
zInv = completeData(:,18:20)*-1;

%% get number of periods of seated rest (events)
numberEvents = completeData(:,2);

%% put data into a table
data = horzcat(subjectIDlist, numClusters', zInv, numberEvents);
columnNames = {'PPID','numClusters','zInv_Neg','zInv_Pos','zInv_M','numEvents'};
data_Table = array2table(data,'VariableNames', columnNames);
save('clusterGranularity_data.mat','data_Table');

%% remove participants with fewer than 7 days of data
index = find(completeData(:,3)<7); % find participants with fewer than minimum number of days
data(index,:) = []; % remove participants with fewer than 7 days of data

%% run correlations
[numClust_mGran_r,numClust_mGran_p] = corr(data(:,2),data(:,5)); % between overall granularity and total number of clusters
[numClust_numEvents_r,numClust_numEvents_p] = corr(data(:,2),data(:,6)); % between total number of clusters and number of seated rests

%% create scatter plot
figure;
scatter1 = scatter(data(:,5),data(:,2),[],rgb('MediumPurple'),'filled');
set(gca,'fontsize',14)
xlim([-1.4 0]);
xlabel('emotional granularity');
ylim([0 14]);
ylabel('number of clusters');
h1 = lsline;
h1.Color = 'k';
saveas(scatter1,'numClust_mGran_scatter_plot','tiff');

%% run multiple regression accounting for number of rests
% model 1: Number of clusters as a function of number of rests and overall granularity
y_numClust = data(:,2); % total number of clusters
m1_X = [ones(size(data,1),1) data(:,6) data(:,5)]; % number of rests + overall granularity
[m1_b,m1_bint,m1_r,m1_rint,m1_stats] = regress(y_numClust,m1_X);
[m1_R2change,m1_Fchange,m1_pchange,m1_df1,m1_df2] = twostep_hierarchical_regression(y_numClust,m1_X,1);

y_numClust_S = zscore(data(:,2)); % z-score variables to get standardized beta as output
m1_X_S = [ones(size(data,1),1) zscore(data(:,6)) zscore(data(:,5))]; % number of rests + overall granularity
[m1_b_S,~,~,~,~] = regress(y_numClust_S,m1_X_S);

% model 0: Number of clusters as a function of number of rests (to get residuals for plotting)
m0_X = [ones(size(data,1),1) data(:,6)]; % number of rests
[~,~,m0_r,~,~] = regress(y_numClust,m0_X);

%% create scatter plot of emotional granularity against residuals
figure;
scatter2 = scatter(data(:,5),m0_r,[],rgb('BlueViolet'),'filled');
set(gca,'fontsize',14)
xlim([-1.4 0]);
xlabel('emotional granularity');
ylabel(['number of clusters (residuals)']);
h1 = lsline;
h1.Color = 'k';
saveas(scatter2,'m1_numClustResid_mGran_scatter_plot','tiff');