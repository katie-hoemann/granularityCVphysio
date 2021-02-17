clear;
clc;

%% load classification and granularity data
load('classificationAnalysis_granularity_data.mat');
N = size(data,1);

%% invert zICC values to get granularity
data(:,2:4) = data(:,2:4)*-1;

%% create data table
data_Table = array2table(data,'VariableNames',{'PPID','zInv_N','zInv_P','zInv_M','emodiv','intens','mPos','mNeg','sdPos','sdNeg','Acc','p'});

%% run correlations
[acc_mGran_r,acc_mGran_p] = corr(data(:,4),data(:,11)); % between overall granularity and classification performance
[acc_pGran_r,acc_pGran_p] = corr(data(:,3),data(:,11)); % between positive granularity and classification performance
[acc_nGran_r,acc_nGran_p] = corr(data(:,2),data(:,11)); % between negative granularity and classification performance

%% create scatter plots
figure;
scatter1 = scatter(data(:,4),data(:,11),[],rgb('SeaGreen'),'filled');
set(gca,'fontsize',14)
xlim([-1.4 0]);
xlabel('(overall) emotional granularity');
%ylim([0 100]);
ylabel('classification accuracy (%)');
h1 = lsline;
h1.Color = 'k';
saveas(scatter1,'acc_mGran_scatter_plot','tiff');

figure;
scatter2 = scatter(data(:,3),data(:,11),[],rgb('DarkGreen'),'filled');
set(gca,'fontsize',14)
xlim([-1.4 0]);
xlabel('positive emotional granularity');
%ylim([0 100]);
ylabel('classification accuracy (%)');
h1 = lsline;
h1.Color = 'k';
saveas(scatter2,'acc_pGran_scatter_plot','tiff');


