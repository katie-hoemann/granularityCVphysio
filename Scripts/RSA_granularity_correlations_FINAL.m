clear all;
clc;

%% load aggregated data
load('seatedRest_results_aggregated_wGranularity.mat'); % load matrix of aggregated event data
headers = completeDataSet.Properties.VariableNames; % get column headers for later reference
data = table2array(completeDataSet); % convert back to array for functions
subjectIDlist = unique(data(:,1)); % get list of unique subject IDs
N = length(subjectIDlist);

%% transform values
data(:,18:20) = data(:,18:20)*-1; % invert zICC values to get granularity (zInv)

%% run zero-order correlations
% between granularity and RSA
[RSA_mGran_r,RSA_mGran_p] = corr(data(:,4),data(:,20),'rows','complete'); % overall granularity
[RSA_pGran_r,RSA_pGran_p] = corr(data(:,4),data(:,19),'rows','complete'); % positive granularity
[RSA_nGran_r,RSA_nGran_p] = corr(data(:,4),data(:,18),'rows','complete'); % negative granularity

%% create scatter plot
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

%% run multiple regressions accounting for RR, IBI, affect
% model 1: RSA as a function of RR and overall granularity
y_RSA = data(:,4); % RSA
m1_X = [ones(size(data,1),1) data(:,16) data(:,20)]; % RR + overall granularity
[m1_b,m1_bint,m1_r,m1_rint,m1_stats] = regress(y_RSA,m1_X);
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

% model 0: RSA as a function of RR (to get residuals for plotting)
m0_X = [ones(size(data,1),1) data(:,16)]; % RR
[~,~,m0_r,~,~] = regress(y_RSA,m0_X);

%% create scatter plot of emotional granularity against residuals
% model 1: granularity against RSA after accounting for RR
figure;
scatter2 = scatter(data(:,20),m0_r,[],rgb('Tomato'),'filled');
set(gca,'fontsize',14)
xlim([-1.4 0]);
xlabel('emotional granularity');
ylim([-2.5 2]);
ylabel('resting RSA (residuals)');
h1 = lsline;
h1.Color = 'k';
saveas(scatter2,'m1_RSAresid_mGran_scatter_plot','tiff');
