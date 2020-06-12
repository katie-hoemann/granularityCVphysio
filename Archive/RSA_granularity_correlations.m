clear all;
clc;

%% load aggregated data (use 'filtered' extension to examine periods where IBI was within 1 SD of participant's mean)
load('seatedRest_results_aggregated_wGranularity.mat'); % load matrix of aggregated event data
headers = completeDataSet.Properties.VariableNames; % get column headers for later reference
data = table2array(completeDataSet); % convert back to array for functions

%% load network measures
load('dataSet18_networkMeasures_Pearson.mat'); % load matrix of network measures per subject

%% load data for day divisions
load('seatedRest_results_dayThirds_RSAvariables.mat'); % thirdsDataCompile

%% set parameters
minDays = 7; % set minimum number of days to be included
numberSDs = 2; % number of standard deviations beyond the mean considered outlier

%% drop participants with only a few days of data
index = find(data(:,3)<minDays); % find participants with fewer than minimum number of days
data(index,:) = []; % remove participants with fewer than minimum number of days
networkMeasures(index,:) = []; % remove participants with fewer than minimum number of days
subjectIDlist = unique(data(:,1)); % get list of unique subject IDs
N = length(subjectIDlist);

%% transform values
data(:,18:20) = data(:,18:20)*-1; % invert zICC values to get granularity (zInv)
data(:,25) = (data(:,4)-mean(data(:,4))).^2; % mean-center and square mean RSA values for testing quadratic fit
data(:,26) = (data(:,14)-mean(data(:,14))).^2; % mean-center and square mean RMSSD values for testing quadratic fit
data(:,27) = (data(:,16)-mean(data(:,16))).^2; % mean-center and square mean RR values for testing quadratic fit
data(:,28) = (data(:,6)-mean(data(:,6))).^2; % mean-center and square mean IBI values for testing quadratic fit
headers = [headers 'RSA_Msq' 'RMSSD_Msq' 'RR_Msq' 'IBI_Msq']; % add new variable names to the header

%% run zero-order correlations
% between granularity and RSA
[RSA_mGran_r,RSA_mGran_p] = corr(data(:,4),data(:,20),'rows','complete'); % overall granularity
[RSA_nGran_r,RSA_nGran_p] = corr(data(:,4),data(:,18),'rows','complete'); % negative granularity
[RSA_pGran_r,RSA_pGran_p] = corr(data(:,4),data(:,19),'rows','complete'); % positive granularity
% between negative diversity and RSA
[RSA_nDiv_r,RSA_nDiv_p] = corr(data(:,4),networkMeasures(:,11),'rows','complete'); % negative diversity coefficient
% between affect and RSA
[RSA_mPos_r,RSA_mPos_p] = corr(data(:,4),data(:,21),'rows','complete'); % mean positive affect
[RSA_sdPos_r,RSA_sdPos_p] = corr(data(:,4),data(:,23),'rows','complete'); % std positive affect
[RSA_mNeg_r,RSA_mNeg_p] = corr(data(:,4),data(:,22),'rows','complete'); % mean negative affect
[RSA_sdNeg_r,RSA_sdNeg_p] = corr(data(:,4),data(:,24),'rows','complete'); % std negative affect

% % broken down by day thirds
% [RSA_mGran_r1,RSA_mGran_p1] = corr(thirdsDataCompile(:,2),data(:,20),'rows','complete'); % overall granularity
% [RSA_mGran_r2,RSA_mGran_p2] = corr(thirdsDataCompile(:,3),data(:,20),'rows','complete'); % overall granularity
% [RSA_mGran_r3,RSA_mGran_p3] = corr(thirdsDataCompile(:,4),data(:,20),'rows','complete'); % overall granularity

% % between granularity and RMSSD
% [RMSSD_mGran_r,RMSSD_mGran_p] = corr(data(:,14),data(:,20),'rows','complete'); % overall granularity
% [RMSSD_nGran_r,RMSSD_nGran_p] = corr(data(:,14),data(:,18),'rows','complete'); % negative granularity
% [RMSSD_pGran_r,RMSSD_pGran_p] = corr(data(:,14),data(:,17),'rows','complete'); % positive granularity
% % between negative diversity and RMSSD
% [RMSSD_nDiv_r,RMSSD_nDiv_p] = corr(data(:,14),networkMeasures(:,11),'rows','complete'); % negative diversity coefficient
% % between affect and RMSSD
% [RMSSD_mPos_r,RMSSD_mPos_p] = corr(data(:,14),data(:,21),'rows','complete'); % mean positive affect
% [RMSSD_sdPos_r,RMSSD_sdPos_p] = corr(data(:,14),data(:,23),'rows','complete'); % std positive affect
% [RMSSD_mNeg_r,RMSSD_mNeg_p] = corr(data(:,14),data(:,22),'rows','complete'); % mean negative affect
% [RMSSD_sdNeg_r,RMSSD_sdNeg_p] = corr(data(:,14),data(:,24),'rows','complete'); % std negative affect

%% create scatter plots to check nature of relationship
figure;
scatter1 = scatter(data(:,20),data(:,4),[],rgb('Coral'),'filled');
set(gca,'fontsize',14)
xlim([-1.4 0]);
xlabel('emotional granularity');
ylim([6 12]);
ylabel('resting RSA');
h1 = lsline;
h1.Color = 'k';
saveas(scatter1,'RSA_mGran_scatter_plot','tiff');

% figure;
% scatter2 = scatter(data(:,20),data(:,14),[],rgb('DarkRed'),'filled');
% xlabel('emotional granularity');
% ylim([25 75]);
% ylabel('resting RMSSD');
% h1 = lsline;
% h1.Color = 'k';
% saveas(scatter2,'RMSSD_mGran_scatter_plot','tiff');

%% check robustness to outliers
m_mGran = mean(data(:,20));
sd_mGran = std(data(:,20));
threshold_mGran_low = m_mGran-(numberSDs*sd_mGran);
threshold_mGran_high = m_mGran+(numberSDs*sd_mGran);
outliers_mGran_low = data(:,20)<threshold_mGran_low;
outliers_mGran_high = data(:,20)>threshold_mGran_high;

m_RSA = mean(data(:,4));
sd_RSA = std(data(:,4));
threshold_RSA_low = m_RSA-(numberSDs*sd_RSA);
threshold_RSA_high = m_RSA+(numberSDs*sd_RSA);
outliers_RSA_low = data(:,4)<threshold_RSA_low;
outliers_RSA_high = data(:,4)>threshold_RSA_high;

outliers_RSA = outliers_mGran_low + outliers_mGran_high + outliers_RSA_low + outliers_RSA_high;
data_noOutliers_RSA = data(~outliers_RSA,:);
numberOutliers_RSA = sum(outliers_RSA);
newN_RSA = size(data_noOutliers_RSA,1);

[RSA_mGran_r_NO,RSA_mGran_p_NO] = corr(data_noOutliers_RSA(:,4),data_noOutliers_RSA(:,20),'rows','complete');

figure;
scatter3 = scatter(data_noOutliers_RSA(:,20),data_noOutliers_RSA(:,4),[],rgb('Blue'),'filled');
xlim([-1.4 0]);
xlabel('emotional granularity');
ylim([6 12]);
ylabel('resting RSA');
h1 = lsline;
h1.Color = 'k';
saveas(scatter3,'RSA_mGran_scatter_plot_noOutliers','tiff');

% m_RMSSD = mean(data(:,14));
% sd_RMSSD = std(data(:,14));
% threshold_RMSSD_low = m_RMSSD-(numberSDs*sd_RMSSD);
% threshold_RMSSD_high = m_RMSSD+(numberSDs*sd_RMSSD);
% outliers_RMSSD_low = data(:,14)<threshold_RMSSD_low;
% outliers_RMSSD_high = data(:,14)>threshold_RMSSD_high;
% 
% outliers_RMSSD = outliers_mGran_low + outliers_mGran_high + outliers_RMSSD_low + outliers_RMSSD_high;
% data_noOutliers_RMSSD = data(~outliers_RMSSD,:);
% numberOutliers_RMSSD = sum(outliers_RMSSD);
% newN_RMSSD = size(data_noOutliers_RMSSD,1);
% 
% [RMSSD_mGran_r_NO,RMSSD_mGran_p_NO] = corr(data_noOutliers_RMSSD(:,14),data_noOutliers_RMSSD(:,20),'rows','complete');
% 
% figure;
% scatter4 = scatter(data_noOutliers_RMSSD(:,20),data_noOutliers_RMSSD(:,14),[],rgb('DarkBlue'),'filled');
% xlim([-1.4 0]);
% xlabel('emotional granularity');
% ylim([25 75]);
% ylabel('resting RMSSD');
% h1 = lsline;
% h1.Color = 'k';
% saveas(scatter2,'RMSSD_mGran_scatter_plot_noOutliers','tiff');

%% run multiple regressions accounting for RR, IBI, affect
% model 1: RSA as a function of RR and overall granularity
y_RSA = data(:,4); % RSA
m1_X = [ones(size(data,1),1) data(:,16) data(:,20)]; % RR + overall granularity
[m1_b,m1_bint,m1_r,m1_rint,m1_stats] = regress(y_RSA,m1_X);
% contain0 = (m1_rint(:,1)<0 & m1_rint(:,2)>0);
% index2 = find(contain0==false);
% y_RSA(index2) = [];
% X_RSA(index2,:) = [];
% [m1_b,m1_bint,m1_r,m1_rint,m1_stats] = regress(y_RSA,m1_X);
[m1_R2change,m1_Fchange,m1_pchange,m1_df1,m1_df2] = twostep_hierarchical_regression(y_RSA,m1_X,1);

y_RSA_S = zscore(data(:,4)); % z-score variables to get standardized beta as output
m1_X_S = [ones(size(data,1),1) zscore(data(:,16)) zscore(data(:,20))]; % RR + overall granularity
[m1_b_S,~,~,~,~] = regress(y_RSA_S,m1_X_S);

% model 2: RSA as a function of RR, IBI, overall granularity
m2_X = [ones(size(data,1),1) data(:,16) data(:,6) data(:,20)]; % RR + IBI + overall granularity
[m2_b,m2_bint,m2_r,m2_rint,m2_stats] = regress(y_RSA,m2_X);
[m2_R2change,m2_Fchange,m2_pchange,m2_df1,m2_df2] = twostep_hierarchical_regression(y_RSA,m2_X,1);

m2_X_S = [ones(size(data,1),1) zscore(data(:,16)) zscore(data(:,6)) zscore(data(:,20))]; % RR + IBI + overall granularity
[m2_b_S,~,~,~,~] = regress(y_RSA_S,m2_X_S); % z-score variables to get standardized beta as output

% model 3: RSA as a function of RR, IBI, positive affect, overall granularity
m3_X = [ones(size(data,1),1) data(:,16) data(:,6) data(:,21) data(:,20)]; % RR + IBI + positive affect + overall granularity
[m3_b,m3_bint,m3_r,m3_rint,m3_stats] = regress(y_RSA,m3_X);
[m3_R2change,m3_Fchange,m3_pchange,m3_df1,m3_df2] = twostep_hierarchical_regression(y_RSA,m3_X,1);

% model 0: RSA as a function of RR (to get residuals for plotting)
m0_X = [ones(size(data,1),1) data(:,16)]; % RR
[~,~,m0_r,~,~] = regress(y_RSA,m0_X);

%% create scatter plots of emotional granularity against residuals
% model 1: granularity against RSA after accounting for RR
figure;
scatter5 = scatter(data(:,20),m0_r,[],rgb('Tomato'),'filled');
set(gca,'fontsize',14)
xlim([-1.4 0]);
xlabel('emotional granularity');
ylim([-2.5 2]);
ylabel('resting RSA (residuals)');
h1 = lsline;
h1.Color = 'k';
saveas(scatter5,'m1_RSAresid_mGran_scatter_plot','tiff');

% model 2: granularity against RSA after accounting for RR, IBI
figure;
scatter6 = scatter(data(:,20),m1_r,[],rgb('DarkTurquoise'),'filled');
xlim([-1.4 0]);
xlabel('emotional granularity');
ylim([-2.5 2]);
ylabel('resting RSA (residuals)');
h1 = lsline;
h1.Color = 'k';
saveas(scatter6,'m2_RSAresid_mGran_scatter_plot','tiff');

% model 3: granularity against RSA after accounting for RR, IBI, positive affect
scatter7 = scatter(data(:,20),m2_r,[],rgb('YellowGreen'),'filled');
xlim([-1.4 0]);
xlabel('emotional granularity');
ylim([-2.5 2]);
ylabel('resting RSA (residuals)');
h1 = lsline;
h1.Color = 'k';
saveas(scatter7,'m2_RSAresid_mGran_scatter_plot','tiff');

%% fit linear mixed-effects models (MLM)
% https://www.mathworks.com/help/stats/fitlme.html
load('seatedRest_results_aggregated_long.mat'); % subjectData_agg
dataTable = array2table(subjectData_agg,'VariableNames',{'PPID','DayID','NumEvents','RSA','RR','IBI','RMSSD','zInv_M','zInv_N','zInv_P'});

% model 1: RSA is DV; RR, mean gran are fixed effects; PP random effects
lme1 = fitlme(dataTable,'RSA~RR+zInv_M+(1|PPID)')
