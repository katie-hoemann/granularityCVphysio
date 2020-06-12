clear all;
clc;

%% load data
load('numEvents_daily.mat'); % 'numEvents_array'
load('RSA_M_daily.mat'); % 'RSA_M_array'
load('RR_M_daily.mat'); % 'RR_M_array'
load('IBI_M_daily.mat'); % 'IBI_M_array'
load('RMSSD_M_daily.mat'); % 'RMSSD_M_array'
load('zInv_M_daily.mat'); % 'zInv_M_array'
load('zInv_Pos_daily.mat'); % 'zInv_Pos_array'
load('zInv_Neg_daily.mat'); % 'zInv_Neg_array'
subjectIDlist = RSA_M_array(2:end,1); % get list of unique subject IDs
N = length(subjectIDlist);

%% set parameters
plotData = 0; % set to 1 to generate scatter plots for each participant
checkMinEvents = 0; % set to 1 to check number of events per day for inclusion
minEvents = 30; % set minimum number of events per day to be included (only applicable if checkMinEvents = 1)
minDays = 7; % set minimum number of days to be included
numberSDs = 2; % number of standard deviations beyond the mean considered outlier
sigThresh = .10; % set threshold for two-tailed significance
trendThresh = .20; % set threshold for trending two-tailed significance

%% run through analyses per participant
for i_subject = 1:length(subjectIDlist)
    %% grab, clean, and prep data
    subjectData = [];
    subjectID = subjectIDlist(i_subject);
    index = find(RSA_M_array(:,1)==subjectID);
    subjectData(:,1) = numEvents_array(index,2:end);
    subjectData(:,2) = RSA_M_array(index,2:end);
    subjectData(:,3) = RR_M_array(index,2:end);
    subjectData(:,4) = IBI_M_array(index,2:end);
    subjectData(:,5) = RMSSD_M_array(index,2:end);
    subjectData(:,6) = zInv_M_array(index,2:end);
    subjectData(:,7) = zInv_Neg_array(index,2:end);
    subjectData(:,8) = zInv_Pos_array(index,2:end);
    % build array for use in MLM analyses
    subjectData_build = [];
    subjectData_build = [subjectID*ones(size(subjectData,1),1) transpose(1:1:size(subjectData,1)) subjectData];
    if i_subject == 1
        subjectData_agg = subjectData_build;
    else
        subjectData_agg = vertcat(subjectData_agg,subjectData_build);
    end
    % remove missing data
    missingData = any(isnan(subjectData),2); 
    subjectData = subjectData(~missingData,:); 
    % remove days without the minimum number of events
    if checkMinEvents == 1
        lowEvents = subjectData(:,1)<minEvents;
        subjectData = subjectData(~lowEvents,:);
    end
    % check to see if subject has minimum number of days
    if size(subjectData,1)<minDays
        continue
    else
        numDays(i_subject) = size(subjectData,1);
    end
    % add mean-centered and squared values for testing quadratic fit
    subjectData(:,9) = (subjectData(:,2)-mean(subjectData(:,2))).^2; % RSA
    subjectData(:,10) = (subjectData(:,3)-mean(subjectData(:,3))).^2; % RR
    subjectData(:,11) = (subjectData(:,4)-mean(subjectData(:,4))).^2; % IBI
    subjectData(:,12) = (subjectData(:,5)-mean(subjectData(:,5))).^2; % RMSSD
    
    %% run zero-order correlations between RSA and granularity
    [RSA_mGran_r(i_subject),RSA_mGran_p(i_subject)] = corr(subjectData(:,2),subjectData(:,6),'rows','complete'); % overall granularity
    [RSA_nGran_r(i_subject),RSA_nGran_p(i_subject)] = corr(subjectData(:,2),subjectData(:,7),'rows','complete'); % negative granularity
    [RSA_pGran_r(i_subject),RSA_pGran_p(i_subject)] = corr(subjectData(:,2),subjectData(:,8),'rows','complete'); % positive granularity
    % create scatter plot to check nature of relationship
    if plotData == 1
        figure;
        scatter1 = scatter(subjectData(:,6),subjectData(:,2),[],rgb('Purple'),'filled');
        xlim([-1.4 0]);
        xlabel('daily emotional granularity');
        ylim([6 12]);
        ylabel('daily resting RSA');
        h1 = lsline;
        h1.Color = 'k';
        title(['PP' num2str(subjectID)]);
        saveas(scatter1,['PP' num2str(subjectID) '_RSA_mGran_scatter_plot'],'tiff');
    end
    
    %% run multiple regressions accounting for RR and IBI
    % model 1: RSA as a function of RR and overall granularity
    y_RSA = subjectData(:,2); % RSA
    m1_X = [ones(size(subjectData,1),1) subjectData(:,3) subjectData(:,6)]; % RR + overall granularity
    [m1_b(:,i_subject),m1_bint(:,:,i_subject),m1_r,~,m1_stats(:,i_subject)] = regress(y_RSA,m1_X);
    [m1_R2change(i_subject),m1_Fchange(i_subject),m1_pchange(i_subject),m1_df1(i_subject),m1_df2(i_subject)] = twostep_hierarchical_regression(y_RSA,m1_X,1);

    % model 2: RSA as a function of RR, IBI, overall granularity
    m2_X = [ones(size(subjectData,1),1) subjectData(:,3) subjectData(:,4) subjectData(:,6)]; % RR + IBI + overall granularity
    [m2_b(:,i_subject),m2_bint(:,:,i_subject),m2_r,~,m2_stats(:,i_subject)] = regress(y_RSA,m2_X);
    [m2_R2change(i_subject),m2_Fchange(i_subject),m2_pchange(i_subject),m2_df1(i_subject),m2_df2(i_subject)] = twostep_hierarchical_regression(y_RSA,m2_X,1);

    % model 0: RSA as a function of RR (to get residuals for plotting)
    m0_X = [ones(size(subjectData,1),1) subjectData(:,3)]; % RR
    [~,~,m0_r,~,~] = regress(y_RSA,m0_X);

    % create scatter plots of emotional granularity against residuals
    if plotData == 1
        figure;
        scatter2 = scatter(subjectData(:,6),m0_r,[],rgb('Orange'),'filled'); % model 1
        xlim([-1.4 0]);
        xlabel('daily emotional granularity');
        ylim([-2.5 2]);
        ylabel('daily resting RSA (residuals)');
        h1 = lsline;
        h1.Color = 'k';
        title(['PP' num2str(subjectID)]);
        saveas(scatter2,['PP' num2str(subjectID) '_m1_RSAresid_mGran_scatter_plot'],'tiff');

        figure;
        scatter3 = scatter(subjectData(:,6),m1_r,[],rgb('Turquoise'),'filled'); % model 2
        xlim([-1.4 0]);
        xlabel('daily emotional granularity');
        ylim([-2.5 2]);
        ylabel('daily resting RSA (residuals)');
        h1 = lsline;
        h1.Color = 'k';
        title(['PP' num2str(subjectID)]);
        saveas(scatter3,['PP' num2str(subjectID) '_m2_RSAresid_mGran_scatter_plot'],'tiff');

        close all;
    end
    clear y_RSA m0_X m1_X m2_X m0_r m1_r m2_r;
end

%% summarize and write results
% zero-order correlations
RSA_mGran = horzcat(subjectIDlist, RSA_mGran_r', RSA_mGran_p', numDays'); % compile correlation results into a single matrix
RSA_nGran = horzcat(subjectIDlist, RSA_nGran_r', RSA_nGran_p', numDays');
RSA_pGran = horzcat(subjectIDlist, RSA_pGran_r', RSA_pGran_p', numDays');
RSA_mGran(RSA_mGran(:,2)==0,:) = []; % remove skipped participants
RSA_nGran(RSA_nGran(:,2)==0,:) = [];
RSA_pGran(RSA_pGran(:,2)==0,:) = [];
newN = length(RSA_mGran); % find final N

RSA_mGran_Table = array2table(RSA_mGran,'VariableNames',{'PPID','r','p','nDays'}); % format table for review/export
RSA_nGran_Table = array2table(RSA_nGran,'VariableNames',{'PPID','r','p','nDays'});
RSA_pGran_Table = array2table(RSA_pGran,'VariableNames',{'PPID','r','p','nDays'});
writetable(RSA_mGran_Table,'RSA_mGran_correlation_summary.xlsx'); % write table
writetable(RSA_nGran_Table,'RSA_nGran_correlation_summary.xlsx');
writetable(RSA_pGran_Table,'RSA_pGran_correlation_summary.xlsx');

RSA_mGran_pos = sum(RSA_mGran(:,2)>0)/newN; % proportion of participants with positive correlations
RSA_mGran_neg = sum(RSA_mGran(:,2)<0)/newN; % proportion of participants with negative correlations
RSA_mGran_sig_pos = sum(RSA_mGran(:,3)<sigThresh & RSA_mGran(:,1)>0)/newN; % proportion of participants with significant positive correlations
RSA_mGran_sig_neg = sum(RSA_mGran(:,3)<sigThresh & RSA_mGran(:,1)<0)/newN; % proportion of participants with significant negative correlations
RSA_mGran_trend_pos = sum(RSA_mGran(:,3)<trendThresh & RSA_mGran(:,1)>0)/newN; % proportion of participants with trending positive correlations
RSA_mGran_trend_neg = sum(RSA_mGran(:,3)<trendThresh & RSA_mGran(:,1)<0)/newN; % proportion of participants with trending negative correlations

RSA_nGran_pos = sum(RSA_nGran(:,2)>0)/newN;
RSA_nGran_neg = sum(RSA_nGran(:,2)<0)/newN;
RSA_nGran_sig_pos = sum(RSA_nGran(:,3)<sigThresh & RSA_nGran(:,1)>0)/newN;
RSA_nGran_sig_neg = sum(RSA_nGran(:,3)<sigThresh & RSA_nGran(:,1)<0)/newN;
RSA_nGran_trend_pos = sum(RSA_nGran(:,3)<trendThresh & RSA_nGran(:,1)>0)/newN;
RSA_nGran_trend_neg = sum(RSA_nGran(:,3)<trendThresh & RSA_nGran(:,1)<0)/newN;

RSA_pGran_pos = sum(RSA_pGran(:,2)>0)/newN;
RSA_pGran_neg = sum(RSA_pGran(:,2)<0)/newN;
RSA_pGran_sig_pos = sum(RSA_pGran(:,3)<sigThresh & RSA_pGran(:,1)>0)/newN;
RSA_pGran_sig_neg = sum(RSA_pGran(:,3)<sigThresh & RSA_pGran(:,1)<0)/newN;
RSA_pGran_trend_pos = sum(RSA_pGran(:,3)<trendThresh & RSA_pGran(:,1)>0)/newN;
RSA_pGran_trend_neg = sum(RSA_pGran(:,3)<trendThresh & RSA_pGran(:,1)<0)/newN;

% hierarchical regressions
m1 = horzcat(subjectIDlist, m1_b(end,:)', m1_R2change', m1_Fchange', m1_pchange', numDays');
m2 = horzcat(subjectIDlist, m2_b(end,:)', m2_R2change', m2_Fchange', m2_pchange', numDays');
m1(m1(:,2)==0,:) = [];
m2(m2(:,2)==0,:) = [];

m1_Table = array2table(m1,'VariableNames',{'PPID','mGran_b','R2change','Fchange','pchange','nDays'});
m2_Table = array2table(m2,'VariableNames',{'PPID','mGran_b','R2change','Fchange','pchange','nDays'});
writetable(m1_Table,'RSA_mGran_model1_regression_summary.xlsx');
writetable(m2_Table,'RSA_mGran_model2_regression_summary.xlsx');

m1_pos = sum(m1(:,2)>0)/newN; % proportion of participants with granularity coefficients that are positive
m1_neg = sum(m1(:,2)<0)/newN; % proportion of participants with granularity coefficients that are negative
m1_sig_pos = sum(m1(:,5)<sigThresh & m1(:,1)>0)/newN; % prop of participants with granularity coefs that are pos and sig improve model fit
m1_sig_neg = sum(m1(:,5)<sigThresh & m1(:,1)<0)/newN; % prop of participants with granularity coefs that are neg and sig improve model fit
m1_trend_pos = sum(m1(:,5)<trendThresh & m1(:,1)>0)/newN; % prop of participants with granularity coefs that are pos and trend improve model fit
m1_trend_neg = sum(m1(:,5)<trendThresh & m1(:,1)<0)/newN; % prop of participants with granularity coefs that are neg and trend improve model fit

m2_pos = sum(m2(:,2)>0)/newN;
m2_neg = sum(m2(:,2)<0)/newN;
m2_sig_pos = sum(m2(:,5)<sigThresh & m2(:,1)>0)/newN;
m2_sig_neg = sum(m2(:,5)<sigThresh & m2(:,1)<0)/newN;
m2_trend_pos = sum(m2(:,5)<trendThresh & m2(:,1)>0)/newN;
m2_trend_neg = sum(m2(:,5)<trendThresh & m2(:,1)<0)/newN;

%% save aggregated long-format data
save('seatedRest_results_aggregated_long.mat','subjectData_agg');
