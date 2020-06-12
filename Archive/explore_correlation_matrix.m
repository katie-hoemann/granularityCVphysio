clear all;
clc;

load('clusterGranularity_data.mat'); % load table of discovered clusters and granularity per participant
data1 = table2array(data_Table); % convert back to array for functions
variableNames1 = data_Table.Properties.VariableNames;
data1_Table = data_Table;

load('seatedRest_results_aggregated_wGranularity.mat'); % load table of aggregated event data
data2 = table2array(completeDataSet); % convert back to array for functions
variableNames2 = completeDataSet.Properties.VariableNames;
data2_Table = completeDataSet;

load('classificationAnalysis_granularity_data.mat'); % load matrix of classification data 
data3 = allData_array;
variableNames3 = {'PPID','zICC_N','zICC_P','zICC_M','emodiv','intens','mPos','mNeg','sdPos','sdNeg','Acc','p'};
data3_Table = array2table(data3,'VariableNames',variableNames3);

load('dataSet18_networkMeasures_Pearson.mat'); % load matrix of network measures
data4 = networkMeasures;
variableNames4 = {'cluster','clusterSD','density','numCom','mod','pPos','pNeg','pPosSD','pNegSD','dPos','dNeg','dPosSD','dNegSD','comRadius','modVA'};
data4_Table = array2table(data4,'VariableNames',variableNames4);

classAcc_Table = outerjoin(data1_Table,data3_Table,'LeftKeys','PPID','RightKeys','PPID');
classAcc = table2array(classAcc_Table(:,18));

varToCorr = [data2(:,[3 13]) data1(:,2) classAcc data1(:,4:6) data4(:,[4 11])];
varToCorr_Names = {'RSA_M','RMSSD_M','numClust','classAcc','zInv_N','zInv_P','zInv_M','numCom','dNeg'};
varToCorr_Table = array2table(varToCorr,'VariableNames',varToCorr_Names);

[varToCorr_r,varToCorr_p] = corr(varToCorr,'rows','pairwise');




