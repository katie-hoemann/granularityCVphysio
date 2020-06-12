clear;
clc;

%% load data
load('classificationAnalysis_granularity_data.mat');
data = allData_array;
N = size(data,1);

%% load number of events submitted to classification
eventData = importdata('classificationAnalysis_events_submitted.xlsx');
eventData.data(eventData.data(:,2)==1,:) = []; % remove participants not used in the analysis
numberEvents = eventData.data(:,3); % get total number events submitted
data = [data numberEvents];

%% set parameters
numberSDs = 2; % number of standard deviations beyond the mean considered outlier

%% invert zICC values to get granularity
data(:,2:4) = data(:,2:4)*-1;
data_Table = array2table(data,'VariableNames',{'PPID','zInv_N','zInv_P','zInv_M','emodiv','intens','mPos','mNeg','sdPos','sdNeg','Acc','p','numEvents'});

%% run correlations
[acc_mGran_r,acc_mGran_p] = corr(data(:,4),data(:,11)); % between overall granularity and classification performance
[acc_pGran_r,acc_pGran_p] = corr(data(:,3),data(:,11)); % between positive granularity and classification performance
[acc_numEvents_r,acc_numEvents_p] = corr(data(:,13),data(:,11)); % between number of events and classification performance

%% create scatter plots
figure;
scatter1 = scatter(data(:,4),data(:,11),[],rgb('SeaGreen'),'filled');
set(gca,'fontsize',14)
xlim([-1.4 0]);
xlabel('emotional granularity');
%ylim([0 100]);
ylabel('classification accuracy (%)');
h1 = lsline;
h1.Color = 'k';
saveas(scatter1,'acc_mGran_scatter_plot','tiff');

% figure;
% scatter2 = scatter(data(:,4),data(:,12),[],rgb('DarkMagenta'),'filled');
% set(gca,'fontsize',14)
% xlabel('emotional granularity');
% ylabel('classification accuracy significance ({\itp})');
% h1 = lsline;
% h1.Color = 'k';
% saveas(scatter2,'sig_mGran_scatter_plot','tiff');

figure;
scatter3 = scatter(data(:,3),data(:,11),[],rgb('DarkGreen'),'filled');
set(gca,'fontsize',14)
xlim([-1.4 0]);
xlabel('positive emotional granularity');
%ylim([0 100]);
ylabel('classification accuracy (%)');
h1 = lsline;
h1.Color = 'k';
saveas(scatter3,'acc_pGran_scatter_plot','tiff');

%% check robustness to outliers
m_mGran = mean(data(:,4));
sd_mGran = std(data(:,4));
threshold_mGran_low = m_mGran-(numberSDs*sd_mGran);
threshold_mGran_high = m_mGran+(numberSDs*sd_mGran);
outliers_mGran_low = data(:,4)<threshold_mGran_low;
outliers_mGran_high = data(:,4)>threshold_mGran_high;

m_acc = mean(data(:,11));
sd_acc = std(data(:,11));
threshold_acc_low = m_acc-(numberSDs*sd_acc);
threshold_acc_high = m_acc+(numberSDs*sd_acc);
outliers_acc_low = data(:,11)<threshold_acc_low;
outliers_acc_high = data(:,11)>threshold_acc_high;

outliers_mGran = outliers_mGran_low + outliers_mGran_high + outliers_acc_low + outliers_acc_high;
data_noOutliers_mGran = data(~outliers_mGran,:);
data_outliers_mGran = data.*outliers_mGran;
data_outliers_mGran(data_outliers_mGran(:,1)==0,:) = [];
numberOutliers_mGran = sum(outliers_mGran);
newN_mGran = size(data_noOutliers_mGran,1);

[acc_mGran_r_NO,acc_mGran_p_NO] = corr(data_noOutliers_mGran(:,4),data_noOutliers_mGran(:,11));
[sig_mGran_r_NO,sig_mGran_p_NO] = corr(data_noOutliers_mGran(:,4),data_noOutliers_mGran(:,12));

figure;
scatter4 = scatter(data_noOutliers_mGran(:,4),data_noOutliers_mGran(:,11),[],rgb('SeaGreen'),'filled');
set(gca,'fontsize',14)
xlim([-1.4 0]);
xlabel('emotional granularity');
ylim([20 90]);
ylabel('classification accuracy (%)');
h1 = lsline;
h1.Color = 'k';
hold on;
scatter(data_outliers_mGran(:,4),data_outliers_mGran(:,11),[],rgb('SeaGreen'));
hold off;
saveas(scatter4,'acc_mGran_scatter_plot_noOutliers','tiff');

m_pGran = mean(data(:,3));
sd_pGran = std(data(:,3));
threshold_pGran_low = m_pGran-(numberSDs*sd_pGran);
threshold_pGran_high = m_pGran+(numberSDs*sd_pGran);
outliers_pGran_low = data(:,3)<threshold_pGran_low;
outliers_pGran_high = data(:,3)>threshold_pGran_high;

outliers_pGran = outliers_pGran_low + outliers_pGran_high + outliers_acc_low + outliers_acc_high;
data_noOutliers_pGran = data(~outliers_pGran,:);
data_outliers_pGran = data.*outliers_pGran;
data_outliers_pGran(data_outliers_pGran(:,1)==0,:) = [];
numberOutliers_pGran = sum(outliers_pGran);
newN_pGran = size(data_noOutliers_pGran,1);

[acc_pGran_r_NO,acc_pGran_p_NO] = corr(data_noOutliers_pGran(:,3),data_noOutliers_pGran(:,11));
[sig_pGran_r_NO,sig_pGran_p_NO] = corr(data_noOutliers_pGran(:,3),data_noOutliers_pGran(:,12));

figure;
scatter5 = scatter(data_noOutliers_pGran(:,3),data_noOutliers_pGran(:,11),[],rgb('SeaGreen'),'filled');
set(gca,'fontsize',14)
xlim([-1.4 0]);
xlabel('positive emotional granularity');
ylim([20 90]);
ylabel('classification accuracy (%)');
h1 = lsline;
h1.Color = 'k';
hold on;
scatter(data_outliers_pGran(:,3),data_outliers_pGran(:,11),[],rgb('SeaGreen'));
hold off;
saveas(scatter5,'acc_pGran_scatter_plot_noOutliers','tiff');

%% run multiple regressions accounting for number of events
y_acc = data(:,11); % classification accuracy
m1_X = [ones(size(data,1),1) data(:,13) data(:,4)]; % number of events + overall granularity
[m1_b,m1_bint,m1_r,m1_rint,m1_stats] = regress(y_acc,m1_X);
[m1_R2change,m1_Fchange,m1_pchange,m1_df1,m1_df2] = twostep_hierarchical_regression(y_acc,m1_X,1);



