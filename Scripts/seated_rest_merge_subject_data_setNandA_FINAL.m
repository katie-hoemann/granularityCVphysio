% created to compile cluster data for DP-GMM sensitivity analyses

clear all;
clc;

%% identify and merge files with results summary per subject
files = dir('C:\Users\Katie\Documents\MATLAB\Seated_rest\Extracted data\**\sitting_features_summary_results_setNandA_*.xlsx'); % specify directory of participant data files to compile
for i_file = 1:length(files) % for every file
    filename = files(i_file).name; % find the name of the file in question
    rawData = importdata(filename); % import the data
    eventData = rawData.data.Sheet1(:,3:29);
    clusterData = rawData.data.Sheet2(:,2:end);
    subjectID = unique(eventData(:,1)); % get subject ID
    subjectIDcluster = subjectID*ones(size(clusterData,1),1); % create subject ID column for cluster data
    clusterData = horzcat(subjectIDcluster,clusterData); % add subject ID column to cluster data
    if i_file == 1
        clusterDataCompile = clusterData; % if it's the first participant, then don't merge with anything
    else
        clusterDataCompile = vertcat(clusterDataCompile,clusterData);
    end
    clear eventData clusterData;
end

%% export data
clusterColumnNames = rawData.textdata.Sheet2(1,:); % grab column headers from cluster data
clusterColumnNames = regexprep(clusterColumnNames,' ','_'); % replace spaces with underscores
clusterColumnNames = ['SID' clusterColumnNames]; % add subjectID column header
clusterDataSet_setNandA = array2table(clusterDataCompile,'VariableNames',clusterColumnNames); % write the cluster data to a table with column headers
writetable(clusterDataSet_setNandA,'cluster_summary_results_setNandA_compiled.xlsx'); % output cluster data to a spreadsheet
save('cluster_results_setNandA_compiled.mat','clusterDataSet_setNandA'); % save the table