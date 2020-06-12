clear all;
clc;

%% set parameters
analysis = 3; % set to the cluster analysis number (e.g., 1 for 'a1')

%% import and clean data
rawData = importdata(['Sitting_features_summary_results_aggregated_clusters_a' num2str(analysis) '.xlsx']);
subjectData = rawData.data.Sheet1(:,2:end);
variableNames = rawData.textdata.Sheet1; % column headers for subject data
clusterData = rawData.data.Sheet2(1:end-1,2:end);
clusterData(:,end) = countcats(categorical(subjectData(:,24))); % get raw number of participants in each cluster
for i_subject = 1:size(subjectData,1) % identify subjects in clusters with fewer than two members
    if any(subjectData(i_subject,24) == find(clusterData(:,end)<2)-1) 
        subjectData(i_subject,24) = NaN;
    end
end
subjectData(isnan(subjectData(:,24))==1,:) = []; % remove subjects excluded from clustering and in small clusters
N = size(subjectData,1);

%% data operations
subjectData(:,17:19) = subjectData(:,17:19)*-1; % multiple zICCs by -1
for i_subject = 1:size(subjectData,1)
    clusterProb(i_subject) = subjectData(i_subject,subjectData(i_subject,24)+25); % get posterior probability of assigned cluster membership
end
subjectData = [subjectData clusterProb'];
subjectData(:,32) = subjectData(:,17).*subjectData(:,31); % weighted negative granularity
subjectData(:,33) = subjectData(:,18).*subjectData(:,31); % weighted positive granularity
subjectData(:,34) = subjectData(:,19).*subjectData(:,31); % weighted mean granularity
subjectData(:,35) = subjectData(:,3).*subjectData(:,31); % weighted mean RSA
subjectData(:,36) = subjectData(:,20).*subjectData(:,31); % weighted mean positive affect
subjectData(:,37) = subjectData(:,21).*subjectData(:,31); % weighted mean negative affect
subjectData(:,38) = subjectData(:,22).*subjectData(:,31); % weighted standard deviation positive affect
subjectData(:,39) = subjectData(:,23).*subjectData(:,31); % weighted standard deviation negative affect
variableNames = [variableNames {'clusterProb','wzInv_N','wzInv_P','wzInv_M','wRSA_M','wmPos','wmNeg','wsdPos','wsdNeg'}]'; % update and transpose header

%% ANOVAs
[p_nGran,tbl_nGran,stats_nGran] = anova1(subjectData(:,32),subjectData(:,24),'off'); % one-way ANOVA testing (weighted) neg gran across clusters
[p_pGran,tbl_pGran,stats_pGran] = anova1(subjectData(:,33),subjectData(:,24),'off'); % one-way ANOVA testing (weighted) pos gran across clusters
[p_mGran,tbl_mGran,stats_mGran] = anova1(subjectData(:,34),subjectData(:,24),'off'); % one-way ANOVA testing (weighted) mean gran across clusters
[p_mRSA,tbl_mRSA,stats_mRSA] = anova1(subjectData(:,35),subjectData(:,24),'off'); % one-way ANOVA testing (weighted) mean RSA across clusters
[p_mPos,tbl_mPos,stats_mPos] = anova1(subjectData(:,36),subjectData(:,24),'off'); % one-way ANOVA testing (weighted) mean pos affect across clusters
[p_mNeg,tbl_mNeg,stats_mNeg] = anova1(subjectData(:,37),subjectData(:,24),'off'); % one-way ANOVA testing (weighted) mean neg affect across clusters
[p_sdPos,tbl_sdPos,stats_sdPos] = anova1(subjectData(:,38),subjectData(:,24),'off'); % one-way ANOVA testing (weighted) SD pos affect across clusters
[p_sdNeg,tbl_sdNeg,stats_sdNeg] = anova1(subjectData(:,39),subjectData(:,24),'off'); % one-way ANOVA testing (weighted) SD neg affect across clusters

%% bar graphs of MI values per analysis
% get data from analysis summaries
rawData_a1 = importdata(['Sitting_features_summary_results_aggregated_clusters_a1.xlsx']);
MIdata_a1 = rawData_a1.data.Sheet2(end,2:end-2);
features_a1 = rawData_a1.textdata.Sheet2(1,2:end-2);
rawData_a2 = importdata(['Sitting_features_summary_results_aggregated_clusters_a2.xlsx']);
MIdata_a2 = rawData_a2.data.Sheet2(end,2:end-2);
features_a2 = rawData_a2.textdata.Sheet2(1,2:end-2);
rawData_a3 = importdata(['Sitting_features_summary_results_aggregated_clusters_a3.xlsx']);
MIdata_a3 = rawData_a3.data.Sheet2(end,2:end-2);
features_a3 = rawData_a3.textdata.Sheet2(1,2:end-2);
rawData_a4 = importdata(['Sitting_features_summary_results_aggregated_clusters_a4.xlsx']);
MIdata_a4 = rawData_a4.data.Sheet2(end,2:end-2);
features_a4 = rawData_a4.textdata.Sheet2(1,2:end-2);

features_a2{5} = 'Gran';
features_a4{9} = 'Gran';
features_a3{2} = 'RSA_{SD}';
features_a3{4} = 'IBI_{SD}';
features_a3{6} = 'RR_{SD}';
features_a3{8} = 'PEP_{SD}';
features_a4{2} = 'RSA_{SD}';
features_a4{4} = 'IBI_{SD}';
features_a4{6} = 'RR_{SD}';
features_a4{8} = 'PEP_{SD}';

features_a1 = categorical(features_a1);
features_a2 = categorical(features_a2);
features_a3 = categorical(features_a3);
features_a4 = categorical(features_a4);

fontsize = 12;
y_ul = .5;
color = rgb('MediumPurple');
%x_label = 'features';
y_label = 'MI';

figure;
% top left bar graph
subplot(2,2,1)
bar_a1 = bar(features_a1,MIdata_a1,'FaceColor',color);
set(gca,'fontsize',fontsize)
ylim([0 y_ul]);
%xlabel(x_label);
ylabel(y_label);
% top right bar graph
subplot(2,2,2)
bar_a2 = bar(features_a2,MIdata_a2,'FaceColor',color);
set(gca,'fontsize',fontsize)
ylim([0 y_ul]);
%xlabel(x_label);
ylabel(y_label);
% bottom left bar graph
subplot(2,2,3)
bar_a3 = bar(features_a3,MIdata_a3,'FaceColor',color);
set(gca,'fontsize',fontsize)
ylim([0 y_ul]);
%xlabel(x_label);
ylabel(y_label);
% bottom right bar graph
subplot(2,2,4)
bar_a4 = bar(features_a4,MIdata_a4,'FaceColor',color);
set(gca,'fontsize',fontsize)
ylim([0 y_ul]);
%xlabel(x_label);
ylabel(y_label);
