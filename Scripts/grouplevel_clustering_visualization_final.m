%% Plot Correlation matrices of feature pairs within clusters
% Specify analysis number at the top, and the rest is automatic
%
% Last modified: Nada K. 05/20/2020

%% 
clear; clc;

analysis_flag = 'a4';
cutoff_numDays = 5;

if strcmp(analysis_flag,'a1')
    measures = {'RSA_M', 'IBI_M', 'RR_M', 'PEP_M'};
    measures_labels = measures;
    subplotsize = [2,3];
elseif strcmp(analysis_flag,'a2')
    measures = {'RSA_M', 'IBI_M', 'RR_M', 'PEP_M','zICC_M'};
    measures_labels = {'RSA_M', 'IBI_M', 'RR_M', 'PEP_M','Gran'};
    subplotsize = [3,3];
elseif strcmp(analysis_flag,'a3')
    measures = {'RSA_M', 'RSA_SD', 'IBI_M', 'IBI_SD','RR_M', 'RR_SD','PEP_M','PEP_SD'};
    measures_labels = measures;
    subplotsize = [2,2];
elseif strcmp(analysis_flag,'a4')
    measures = {'RSA_M', 'RSA_SD', 'IBI_M', 'IBI_SD','RR_M', 'RR_SD','PEP_M','PEP_SD','zICC_M'};
    measures_labels = {'RSA_M', 'RSA_SD', 'IBI_M', 'IBI_SD','RR_M', 'RR_SD','PEP_M','PEP_SD','Gran'};
    subplotsize = [2,3];
else
    error('Oops, try again');
end

% upload spreadsheet:
sheet = readtable(['Sitting_features_summary_results_aggregated_clusters_' analysis_flag '.xlsx']);
sheet2 = readtable(['Sitting_features_summary_results_aggregated_clusters_' analysis_flag '.xlsx'], ...
    'Sheet','Sheet2');
sheet(sheet.numDays<=cutoff_numDays,:) = [];

% subjects data:
subjects_physio = [];
for i = 1:length(measures)
    subjects_physio = [subjects_physio, sheet.(measures{i})];
end

if strcmp(analysis_flag,'a2') || strcmp(analysis_flag,'a4')
    subjects_physio(:,5) = -1.*subjects_physio(:,5); % correct for granularity
end

numClusters = length(unique(sheet.C_ID));
for i =1:numClusters
    clusters_labels{i} = ['Cluster - ' num2str(i-1)];
    clusters_length(i) = sum(sheet.C_ID==(i-1));
end

% Correlation for every cluster 
figure; count = 1;
for i = 1:numClusters 
    if clusters_length(i)==1 % skips clusters with one subject
        continue;
    end
    subplot(subplotsize(1),subplotsize(2),count);
    data2plot = corrcoef(subjects_physio(sheet.C_ID==(i-1),:));
    imagesc(data2plot); title([clusters_labels{i}],'interpreter','none');
    colorbar; caxis([-1 1]); set(gca,'fontsize',16); axis('square');
    ax = gca; ax.XTick = 1:length(measures_labels); ax.YTick = 1:length(measures_labels);
    ax.XTickLabel = measures_labels; ax.YTickLabel = measures_labels; 
    ax.XTickLabelRotation = 45;
    count = count + 1;
end
