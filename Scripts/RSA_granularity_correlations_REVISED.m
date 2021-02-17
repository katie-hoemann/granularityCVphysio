% revised to incorporate outlier detection and RSA adjustment for IBI

clear all;
clc;

%% set parameters
dropOutliers = 1; % set to 1 to drop observations with Cook's D > 4/n

%% load aggregated data
load('seatedRest_results_aggregated_wGranularity.mat'); % load matrix of aggregated event data
headers = completeDataSet.Properties.VariableNames; % get column headers for later reference
data = table2array(completeDataSet); % convert back to array for functions
subjectIDlist = unique(data(:,1)); % get list of unique subject IDs

%% pull out specific variables, transform values
granularity = data(:,19:21)*-1; % invert zICC values to get granularity (zInv)
RSA = data(:,4);
RR = data(:,16);
IBI = data(:,6);
RSA_adjusted = 100*(RSA./log(IBI.*IBI)); % adjust RSA according to the formula from de Geus et al (2019): 100*[ln(HFpower)/ln(IBI*IBI)]; ln(HFpower) = RSA

%% run zero-order correlations
% between granularity and RSA
[RSA_mGran_r,RSA_mGran_p] = corr(RSA,granularity(:,3),'rows','complete'); % overall granularity

% between granularity and IBI
[IBI_mGran_r,IBI_mGran_p] = corr(IBI,granularity(:,3),'rows','complete'); % overall granularity

%% run multiple regressions controlling for RR, adjusting for IBI
% model 1: RSA as a function of RR and overall granularity
y1 = RSA; 
predictors = [RR granularity(:,3)]; 
mdl1 = fitlm(predictors,y1);
if dropOutliers == 1
    outliers1 = find(mdl1.Diagnostics.CooksDistance>4/mdl1.NumObservations); % use Cook's distance to identify outliers
    mdl1 = fitlm(predictors,y1,'exclude',outliers1); % re-run model excluding outliers
end

% z-score variables to get standardized beta as output
y1_S = zscore(RSA); 
predictors_S = [zscore(RR) zscore(granularity(:,3))]; 
mdl1_S = fitlm(predictors_S,y1_S);
if dropOutliers == 1
    mdl1_S = fitlm(predictors_S,y1_S,'exclude',outliers1);
end

% use alternative multiple regression functions to get R^2 and F stats
if dropOutliers == 1
    y1_NO = y1;
    y1_NO(outliers1) = [];
    RR_NO = RR;
    RR_NO(outliers1) = [];
    granularity_NO = granularity;
    granularity_NO(outliers1,:) = [];
    predictors_padded_NO = [ones(length(RR_NO),1) RR_NO granularity_NO(:,3)]; 
    [mdl1_b,mdl1_bint,~,~,~] = regress(y1_NO,predictors_padded_NO);
    [mdl1_R2change,mdl1_Fchange,mdl1_pchange,mdl1_df1,mdl1_df2] = twostep_hierarchical_regression(y1_NO,predictors_padded_NO,1);
else
    predictors_padded = [ones(length(RR),1) RR granularity(:,3)]; 
    [mdl1_b,mdl1_bint,~,~,~] = regress(y1,predictors_padded);
    [mdl1_R2change,mdl1_Fchange,mdl1_pchange,mdl1_df1,mdl1_df2] = twostep_hierarchical_regression(y1,predictors_padded,1);
end

% model RSA as a function of RR (to get residuals for plotting)
predictors_R = RR; 
mdl_R = fitlm(predictors_R,y1);
y1_residuals = table2array(mdl_R.Residuals(:,1));
if dropOutliers == 1
    y1_residuals_NO = y1_residuals;
    y1_residuals_NO(outliers1,:) = [];
end

% model 2: adjusted RSA as a function of RR and overall granularity
y2 = RSA_adjusted;
mdl2 = fitlm(predictors,y2);
if dropOutliers == 1
    outliers2 = find(mdl2.Diagnostics.CooksDistance>4/mdl2.NumObservations);
    mdl2 = fitlm(predictors,y2,'exclude',outliers2);
end

% z-score variables to get standardized beta as output
y2_S = zscore(RSA_adjusted); 
mdl2_S = fitlm(predictors_S,y2_S);
if dropOutliers == 1
    mdl2_S = fitlm(predictors_S,y2_S,'exclude',outliers2);
end

% use alternative multiple regression functions to get R^2 and F stats
if dropOutliers == 1
    y2_NO = y2;
    y2_NO(outliers1) = [];
    [mdl2_b,mdl2_bint,~,~,~] = regress(y2_NO,predictors_padded_NO);
    [mdl2_R2change,mdl2_Fchange,mdl2_pchange,mdl2_df1,mdl2_df2] = twostep_hierarchical_regression(y2_NO,predictors_padded_NO,1);
else
    [mdl2_b,mdl2_bint,~,~,~] = regress(y2,predictors_padded);
    [mdl2_R2change,mdl2_Fchange,mdl2_pchange,mdl2_df1,mdl2_df2] = twostep_hierarchical_regression(y2,predictors_padded,1);
end

%% create summary tables
mdl1_summary = [mdl1.NumObservations; mdl1_b(end); mdl1_bint(end,1); mdl1_bint(end,2); mdl1_S.Coefficients.Estimate(end);...
    mdl1.Coefficients.tStat(end); mdl1.DFE; mdl1_R2change; mdl1_Fchange; mdl1_df1; mdl1_df2; mdl1_pchange/2];
mdl1_summary_Table = array2table(mdl1_summary,'RowNames',{'N','b','b_CI_l','b_CI_u','B','t','t_df','R2','F','F_df1','F_df2','p'});

mdl2_summary = [mdl2.NumObservations; mdl2_b(end); mdl2_bint(end,1); mdl2_bint(end,2); mdl2_S.Coefficients.Estimate(end);...
    mdl2.Coefficients.tStat(end); mdl2.DFE; mdl2_R2change; mdl2_Fchange; mdl2_df1; mdl2_df2; mdl2_pchange/2];
mdl2_summary_Table = array2table(mdl2_summary,'RowNames',{'N','b','b_CI_l','b_CI_u','B','t','t_df','R2','F','F_df1','F_df2','p'});

%% create scatter plot of emotional granularity against residuals
if dropOutliers == 1
    figure;
    granularity_NO = granularity;
    granularity_NO(outliers1,:) = [];
    y1_residuals_NO = y1_residuals;
    y1_residuals_NO(outliers1,:) = [];
    scatter1 = scatter(granularity_NO(:,3),y1_residuals_NO,[],rgb('CornflowerBlue'),'filled');
    set(gca,'fontsize',14)
    xlim([-1.4 0]);
    xlabel('emotional granularity');
    ylim([-2.5 2]);
    ylabel('resting RSA (residuals)');
    h1 = lsline;
    h1.Color = 'k';
    hold on;
    scatter2 = scatter(granularity(outliers1,3),y1_residuals(outliers1),[],rgb('Tomato'),'filled');
    hold off;
else
    figure;
    scatter1 = scatter(granularity(:,3),y1_residuals,[],rgb('CornflowerBlue'),'filled');
    set(gca,'fontsize',14)
    xlim([-1.4 0]);
    xlabel('emotional granularity');
    ylim([-2.5 2]);
    ylabel('resting RSA (residuals)');
    h1 = lsline;
    h1.Color = 'k';
end