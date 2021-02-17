clear;
clc;

%% set parameters
dropOutliers = 1; % set to 1 to drop observations with Cook's D > 4/n

%% load classification and granularity data
load('classificationAnalysis_granularity_data.mat');

%% pull out variables
granularity = data(:,19:21)*-1; % invert zICC values to get granularity
classAcc = data(:,30); % classification accuracy for linear SVM on 3 features for 3 classes (i.e., emotion words)
classAcc_chance = data(:,31); % classification accuracy using dummy classifier
numEvents = data(:,29); % number of events that could be sampled for all 3 classes

%% examine model performance against chance
[~,t_p,~,t_stats] = ttest(classAcc,classAcc_chance);

%% run zero-order correlations
[acc_mGran_r,acc_mGran_p] = corr(granularity(:,3),classAcc,'rows','complete'); % between overall granularity and classification performance

%% run regression
% model classification performance as a function of overall granularity
y = classAcc; 
predictors = granularity(:,3); 
mdl = fitlm(predictors,y);
if dropOutliers == 1
    outliers = find(mdl.Diagnostics.CooksDistance>4/mdl.NumObservations);
    mdl = fitlm(predictors,y,'exclude',outliers);
end
    
% z-score variables to get standardized beta as output
y_S = (classAcc-nanmean(classAcc))./nanstd(classAcc);
predictors_S = zscore(granularity(:,3)); 
mdl_S = fitlm(predictors_S,y_S);
if dropOutliers == 1
    mdl_S = fitlm(predictors_S,y_S,'exclude',outliers);
end
   
% use alternative multiple regression functions to get R^2 and F stats
if dropOutliers == 1
    y_NO = y;
    y_NO(outliers) = [];
    granularity_NO = granularity;
    granularity_NO(outliers,:) = [];
    predictors_padded_NO = [ones(size(granularity_NO,1),1) granularity_NO(:,3)]; 
    [mdl_b,mdl_bint,~,~,mdl_stats] = regress(y_NO,predictors_padded_NO);
else
    predictors_padded = [ones(size(granularity,1),1) granularity(:,3)]; 
    [mdl_b,mdl_bint,~,~,mdl_stats] = regress(y,predictors_padded);
end

%% create summary table
mdl_summary = [mdl.NumObservations; mdl_b(end); mdl_bint(end,1); mdl_bint(end,2); mdl_S.Coefficients.Estimate(end);...
    mdl.Coefficients.tStat(end); mdl.DFE; mdl_stats(1); mdl_stats(2); mdl.NumPredictors; mdl.DFE; mdl_stats(3)/2];
mdl_summary_Table = array2table(mdl_summary,'RowNames',{'N','b','b_CI_l','b_CI_u','B','t','t_df','R2','F','F_df1','F_df2','p'});

%% create scatter plots
if dropOutliers == 1
    figure;
    scatter1 = scatter(granularity_NO(:,3),y_NO,[],rgb('SeaGreen'),'filled');
    set(gca,'fontsize',14)
    xlim([-1.4 0]);
    xlabel('emotional granularity');
    ylim([.2 .5]);
    ylabel('classification accuracy');
    h1 = lsline;
    h1.Color = 'k';
    hold on;
    scatter2 = scatter(granularity(outliers,3),y(outliers),[],rgb('Tomato'),'filled');
    hold off;
else
    figure;
    scatter1 = scatter(granularity(:,3),y,[],rgb('SeaGreen'),'filled');
    set(gca,'fontsize',14)
    xlim([-1.4 0]);
    xlabel('emotional granularity');
    ylim([.2 .5]);
    ylabel('classification accuracy');
    h1 = lsline;
    h1.Color = 'k';
end