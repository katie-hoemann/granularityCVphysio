clear;
clc;

load('seatedRest_results_aggregated_wGranularity.mat'); % completeDataSet
rawData = importdata('Supervised clustering results_updated.xlsx'); % import updated classification results
rawData_Table = array2table(rawData.data,'VariableNames',rawData.textdata(1,:)); % transform to table for joining
fullData_Table = outerjoin(completeDataSet,rawData_Table,'LeftKeys','PPID','RightKeys','PPID'); % join tables based on PPID
fullData_Table = removevars(fullData_Table,'PPID_rawData_Table'); % remove duplicate variable
fullData_Table.Properties.VariableNames{1} = 'PPID'; % relabel variable
fullData_Table = removevars(fullData_Table,{'Word1','Word2','Word3'}); % remove text variables
data = table2array(fullData_Table); % create matrix for saving data
writetable(fullData_Table,'classificationAnalysis_granularity_data.xlsx'); % output data to a spreadsheet
save('classificationAnalysis_granularity_data.mat','data'); % save the table