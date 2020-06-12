clear all;
clc;

%% set parameters
filterIBI = 0; % set to 1 to filter from statistics events where IBI is more than 1 SD away from participant's grand mean
divideDay = 1; % set to 1 to divide day into 3 groups of events 

%% identify and merge files with results summary per subject
files = dir('C:\Users\Katie\Documents\MATLAB\Seated_rest\Extracted data\**\sitting_features_summary_results_*.xlsx'); % specify directory of participant data files to compile
for i_file = 1:length(files) % for every file
    filename = files(i_file).name; % find the name of the file in question
    rawData = importdata(filename); % import the data
    eventData = rawData.data.Sheet1(:,3:29);
    clusterData = rawData.data.Sheet2(:,2:end);
    subjectID = unique(eventData(:,1)); % get subject ID
    subjectIDcol = subjectID*ones(size(clusterData,1),1); % create subject ID column
    clusterData = horzcat(subjectIDcol,clusterData); % add subject ID column to cluster data
    if i_file == 1
        eventDataCompile = eventData; % if it's the first participant, then don't merge with anything
        clusterDataCompile = clusterData;
    else
        eventDataCompile = vertcat(eventDataCompile,eventData);
        clusterDataCompile = vertcat(clusterDataCompile,clusterData);
    end
    clear eventData clusterData;
end

%% check data, add headers, and re-export
subjectIDlist = unique(eventDataCompile(:,1)); % should be appropriate length vector

eventColumnNames = rawData.textdata.Sheet1(1,2:28); % grab column headers from event data
eventColumnNames = regexprep(eventColumnNames,' ','_'); % replace spaces with underscores
eventColumnNames = regexprep(eventColumnNames,'(.*',''); % remove everything in parentheses
eventDataSet = array2table(eventDataCompile,'VariableNames',eventColumnNames); % write the event data to a table with column headers
writetable(eventDataSet,'sitting_features_summary_results_compiled.xlsx'); % output event data to a spreadsheet
save('seatedRest_results_compiled.mat','eventDataSet'); % save the table

clusterColumnNames = rawData.textdata.Sheet2(1,:); % grab column headers from cluster data
clusterColumnNames = regexprep(clusterColumnNames,' ','_'); % replace spaces with underscores
clusterColumnNames = ['SID' clusterColumnNames]; % add subjectID column header
clusterDataSet = array2table(clusterDataCompile,'VariableNames',clusterColumnNames); % write the cluster data to a table with column headers
writetable(clusterDataSet,'cluster_summary_results_compiled.xlsx'); % output cluster data to a spreadsheet
save('cluster_results_compiled.mat','clusterDataSet'); % save the table

%% filter data
eventDataCompile_clean = eventDataCompile(~eventDataCompile(:,12)==1,:); % drop events flagged as bad

if filterIBI == 1 % drop events where IBI is more than 1 SD away from mean
    for i_subject = 1:length(subjectIDlist)
        subjectData = [];
        subjectID = subjectIDlist(i_subject);
        index = find(eventDataCompile_clean(:,1)==subjectID);
        subjectData = eventDataCompile_clean(index,:);
        IBI_grandMean = mean(subjectData(:,16));
        IBI_grandSD = std(subjectData(:,16));
        flagIBI_subject = zeros(size(subjectData,1),1);
        for i_event = 1:size(subjectData,1);
            if subjectData(i_event,16) > IBI_grandMean+IBI_grandSD
                flagIBI_subject(i_event) = 1;
            elseif subjectData(i_event,16) < IBI_grandMean-IBI_grandSD
                flagIBI_subject(i_event) = 1;
            end
        end
        if i_subject == 1
            flagIBI = flagIBI_subject;
        else
            flagIBI = [flagIBI; flagIBI_subject];
        end
        clear flagIBI_subject
    end
eventDataCompile_clean = eventDataCompile_clean(flagIBI==0,:); 
end

%% get summary statistics per subject and per day
for i_subject = 1:length(subjectIDlist)
    subjectData = [];
    subjectID = subjectIDlist(i_subject);
    index1 = find(eventDataCompile_clean(:,1)==subjectID);
    subjectData = eventDataCompile_clean(index1,:);
    numberEvents(i_subject) = size(subjectData,1);
    dayIDlist = unique(subjectData(:,2));
    numberDays(i_subject) = length(dayIDlist);
    RSA_grandMean(i_subject) = mean(subjectData(:,14));
    RSA_grandSD(i_subject) = std(subjectData(:,14));
    IBI_grandMean(i_subject) = mean(subjectData(:,16));
    IBI_grandSD(i_subject) = std(subjectData(:,16));
    PEP_grandMean(i_subject) = mean(subjectData(:,18));
    PEP_grandSD(i_subject) = std(subjectData(:,18));
    LVET_grandMean(i_subject) = mean(subjectData(:,20));
    LVET_grandSD(i_subject) = std(subjectData(:,20));
    CO_grandMean(i_subject) = mean(subjectData(:,22));
    CO_grandSD(i_subject) = std(subjectData(:,22));
    RMSSD_grandMean(i_subject) = mean(subjectData(:,24));
    RMSSD_grandSD(i_subject) = std(subjectData(:,24));
    RR_grandMean(i_subject) = mean(subjectData(:,26));
    RR_grandSD(i_subject) = std(subjectData(:,26));
    for i_day = 1:length(dayIDlist)
        dayData = [];
        dayID = dayIDlist(i_day);
        index2 = find(subjectData(:,2)==dayID);
        dayData = subjectData(index2,:);
        dayData = sortrows(dayData,3,'ascend');
        numberEvents_day(i_day) = size(dayData,1);
        RSA_dayMean(i_day) = mean(dayData(:,14));
        RSA_daySD(i_day) = std(dayData(:,14));
        IBI_dayMean(i_day) = mean(dayData(:,16));
        IBI_daySD(i_day) = std(dayData(:,16));
        PEP_dayMean(i_day) = mean(dayData(:,18));
        PEP_daySD(i_day) = std(dayData(:,18));
        LVET_dayMean(i_day) = mean(dayData(:,20));
        LVET_daySD(i_day) = std(dayData(:,20));
        CO_dayMean(i_day) = mean(dayData(:,22));
        CO_daySD(i_day) = std(dayData(:,22));
        RMSSD_dayMean(i_day) = mean(dayData(:,24));
        RMSSD_daySD(i_day) = std(dayData(:,24));
        RR_dayMean(i_day) = mean(dayData(:,26));
        RR_daySD(i_day) = std(dayData(:,26));
        if divideDay == 1
            if numberEvents_day(i_day) < 3
                continue
            end
            if mod(numberEvents_day(i_day),3) > 0
                numberEvents_day_rounded = numberEvents_day(i_day)+3-mod(numberEvents_day(i_day),3);
                sizeThirds = numberEvents_day_rounded/3;
            else
                sizeThirds = numberEvents_day(i_day)/3;
            end
            firstThird = dayData(1:sizeThirds,:);
                RSA_thirds(i_day,1) = mean(firstThird(:,14));
                IBI_thirds(i_day,1) = mean(firstThird(:,16));
                RR_thirds(i_day,1) = mean(firstThird(:,26));
            secondThird = dayData(1+sizeThirds:2*sizeThirds,:);
                RSA_thirds(i_day,2) = mean(secondThird(:,14));
                IBI_thirds(i_day,2) = mean(secondThird(:,16));
                RR_thirds(i_day,2) = mean(secondThird(:,26));
            thirdThird = dayData(1+2*sizeThirds:end,:);
                RSA_thirds(i_day,3) = mean(thirdThird(:,14));
                IBI_thirds(i_day,3) = mean(thirdThird(:,16));
                RR_thirds(i_day,3) = mean(thirdThird(:,26));
        end
        clear firstThird secondThird thirdThird
    end
    % aggregate data for day divisions, if applicable
    if divideDay == 1
        RSA_thirdsMean(i_subject,:) = mean(RSA_thirds,1,'omitnan');
        IBI_thirdsMean(i_subject,:) = mean(IBI_thirds,1,'omitnan');
        RR_thirdsMean(i_subject,:) = mean(RR_thirds,1,'omitnan');
    end
    % add data to tables
    ppID = ['PP' num2str(subjectID)];
    if i_subject == 1
        day = (1:1:max(eventDataCompile_clean(:,2)))';
        day = array2table(day,'VariableNames',{'Day'});
                
        numEvents_day = tablecompile_start(day,dayIDlist,numberEvents_day',ppID,'Day','ppDay','ppDay');
        RSA_M = tablecompile_start(day,dayIDlist,RSA_dayMean',ppID,'Day','ppDay','ppDay');
        RSA_SD = tablecompile_start(day,dayIDlist,RSA_daySD',ppID,'Day','ppDay','ppDay');
        IBI_M = tablecompile_start(day,dayIDlist,IBI_dayMean',ppID,'Day','ppDay','ppDay');
        IBI_SD = tablecompile_start(day,dayIDlist,IBI_daySD',ppID,'Day','ppDay','ppDay');
        PEP_M = tablecompile_start(day,dayIDlist,PEP_dayMean',ppID,'Day','ppDay','ppDay');
        PEP_SD = tablecompile_start(day,dayIDlist,PEP_daySD',ppID,'Day','ppDay','ppDay');
        LVET_M = tablecompile_start(day,dayIDlist,LVET_dayMean',ppID,'Day','ppDay','ppDay');
        LVET_SD = tablecompile_start(day,dayIDlist,LVET_daySD',ppID,'Day','ppDay','ppDay');
        CO_M = tablecompile_start(day,dayIDlist,CO_dayMean',ppID,'Day','ppDay','ppDay');
        CO_SD = tablecompile_start(day,dayIDlist,CO_daySD',ppID,'Day','ppDay','ppDay');
        RMSSD_M = tablecompile_start(day,dayIDlist,RMSSD_dayMean',ppID,'Day','ppDay','ppDay');
        RMSSD_SD = tablecompile_start(day,dayIDlist,RMSSD_daySD',ppID,'Day','ppDay','ppDay');
        RR_M = tablecompile_start(day,dayIDlist,RR_dayMean',ppID,'Day','ppDay','ppDay');
        RR_SD = tablecompile_start(day,dayIDlist,RR_daySD',ppID,'Day','ppDay','ppDay');
    else
        numEvents_day = tablecompile_iter(numEvents_day,dayIDlist,numberEvents_day',ppID,'Day','ppDay','ppDay');
        RSA_M = tablecompile_iter(RSA_M,dayIDlist,RSA_dayMean',ppID,'Day','ppDay','ppDay');
        RSA_SD = tablecompile_iter(RSA_SD,dayIDlist,RSA_daySD',ppID,'Day','ppDay','ppDay');
        IBI_M = tablecompile_iter(IBI_M,dayIDlist,IBI_dayMean',ppID,'Day','ppDay','ppDay');
        IBI_SD = tablecompile_iter(IBI_SD,dayIDlist,IBI_daySD',ppID,'Day','ppDay','ppDay');
        PEP_M = tablecompile_iter(PEP_M,dayIDlist,PEP_dayMean',ppID,'Day','ppDay','ppDay');
        PEP_SD = tablecompile_iter(PEP_SD,dayIDlist,PEP_daySD',ppID,'Day','ppDay','ppDay');
        LVET_M = tablecompile_iter(LVET_M,dayIDlist,LVET_dayMean',ppID,'Day','ppDay','ppDay');
        LVET_SD = tablecompile_iter(LVET_SD,dayIDlist,LVET_daySD',ppID,'Day','ppDay','ppDay');
        CO_M = tablecompile_iter(CO_M,dayIDlist,CO_dayMean',ppID,'Day','ppDay','ppDay');
        CO_SD = tablecompile_iter(CO_SD,dayIDlist,CO_daySD',ppID,'Day','ppDay','ppDay');
        RMSSD_M = tablecompile_iter(RMSSD_M,dayIDlist,RMSSD_dayMean',ppID,'Day','ppDay','ppDay');
        RMSSD_SD = tablecompile_iter(RMSSD_SD,dayIDlist,RMSSD_daySD',ppID,'Day','ppDay','ppDay');
        RR_M = tablecompile_iter(RR_M,dayIDlist,RR_dayMean',ppID,'Day','ppDay','ppDay');
        RR_SD = tablecompile_iter(RR_SD,dayIDlist,RR_daySD',ppID,'Day','ppDay','ppDay');
    end
    clear numberEvents_day RSA_dayMean RSA_daySD IBI_dayMean IBI_daySD PEP_dayMean PEP_daySD LVET_dayMean LVET_daySD CO_dayMean CO_daySD RMSSD_dayMean RMSSD_daySD RR_dayMean RR_daySD;
end

%% compile, write, and save data
% write and save subject-level data
eventDataCompile_agg = [subjectIDlist numberEvents' numberDays' RSA_grandMean' RSA_grandSD' IBI_grandMean' IBI_grandSD' PEP_grandMean' PEP_grandSD'...
    LVET_grandMean' LVET_grandSD' CO_grandMean' CO_grandSD' RMSSD_grandMean' RMSSD_grandSD' RR_grandMean' RR_grandSD'];
aggColumnNames = {'PPID' 'numEvents' 'numDays' 'RSA_M' 'RSA_SD' 'IBI_M' 'IBI_SD' 'PEP_M' 'PEP_SD' 'LVET_M' 'LVET_SD' 'CO_M' 'CO_SD' 'RMSSD_M' 'RMSSD_SD' 'RR_M' 'RR_SD'};
aggDataSet = array2table(eventDataCompile_agg,'VariableNames',aggColumnNames); % write the data to a table with column headers
if filterIBI == 1
    writetable(aggDataSet,'sitting_features_summary_results_aggregated_filtered.xlsx'); % output data to a spreadsheet
    save('seatedRest_results_aggregated_filtered.mat','aggDataSet'); % save the table
else
    writetable(aggDataSet,'sitting_features_summary_results_aggregated.xlsx'); % output data to a spreadsheet
    save('seatedRest_results_aggregated.mat','aggDataSet'); % save the table
end

% add granularity variable(s) of interest
load('granularity_affect_measures.mat'); % load matrix of granularity variables per subject
variableNames = {'PPID' 'zICC_N','zICC_P','zICC_M','mPos','mNeg','sdPos','sdNeg'}; 
subjectMeasures_Table = array2table(subjectMeasures,'VariableNames',variableNames); % convert to table for joining to data set
completeDataSet = outerjoin(aggDataSet,subjectMeasures_Table,'LeftKeys','PPID','RightKeys','PPID'); % join tables based on PPID
completeDataSet = removevars(completeDataSet,'PPID_subjectMeasures_Table'); % remove duplicate variable
completeDataSet.Properties.VariableNames{1} = 'PPID'; % relabel variable
if filterIBI == 1
    writetable(completeDataSet,'sitting_features_summary_results_aggregated_w granularity_filtered.xlsx'); % output data to a spreadsheet
    save('seatedRest_results_aggregated_wGranularity_filtered.mat','completeDataSet'); % save the table
else
    writetable(completeDataSet,'sitting_features_summary_results_aggregated_w granularity.xlsx'); % output data to a spreadsheet
    save('seatedRest_results_aggregated_wGranularity.mat','completeDataSet'); % save the table
end

% save data for day divisions, if applicable
if divideDay == 1
    thirdsDataCompile = [subjectIDlist RSA_thirdsMean IBI_thirdsMean RR_thirdsMean];
    save('seatedRest_results_dayThirds_RSAvariables.mat','thirdsDataCompile'); % save the array
end

% write day-level data
writetable(numEvents_day,'numEvents_daily.xlsx');
writetable(RSA_M,'RSA_M_daily.xlsx');
writetable(RSA_SD,'RSA_SD_daily.xlsx');
writetable(IBI_M,'IBI_M_daily.xlsx');
writetable(IBI_SD,'IBI_SD_daily.xlsx');
writetable(PEP_M,'PEP_M_daily.xlsx');
writetable(PEP_SD,'PEP_SD_daily.xlsx');
writetable(LVET_M,'LVET_M_daily.xlsx');
writetable(LVET_SD,'LVET_SD_daily.xlsx');
writetable(CO_M,'CO_M_daily.xlsx');
writetable(CO_SD,'CO_SD_daily.xlsx');
writetable(RMSSD_M,'RMSSD_M_daily.xlsx');
writetable(RMSSD_SD,'RMSSD_SD_daily.xlsx');
writetable(RR_M,'RR_M_daily.xlsx');
writetable(RR_SD,'RR_SD_daily.xlsx');

% save day-level data
subjectCol = [0; subjectIDlist];

numEvents_array = table2array(numEvents_day)';
numEvents_array = horzcat(subjectCol,numEvents_array);
save('numEvents_daily.mat','numEvents_array');

RSA_M_array = table2array(RSA_M)';
RSA_M_array = horzcat(subjectCol,RSA_M_array);
save('RSA_M_daily.mat','RSA_M_array');

RSA_SD_array = table2array(RSA_SD)';
RSA_SD_array = horzcat(subjectCol,RSA_SD_array);
save('RSA_SD_daily.mat','RSA_SD_array');

IBI_M_array = table2array(IBI_M)';
IBI_M_array = horzcat(subjectCol,IBI_M_array);
save('IBI_M_daily.mat','IBI_M_array');

IBI_SD_array = table2array(IBI_SD)';
IBI_SD_array = horzcat(subjectCol,IBI_SD_array);
save('IBI_SD_daily.mat','IBI_SD_array');

PEP_M_array = table2array(PEP_M)';
PEP_M_array = horzcat(subjectCol,PEP_M_array);
save('PEP_M_daily.mat','PEP_M_array');

PEP_SD_array = table2array(PEP_SD)';
PEP_SD_array = horzcat(subjectCol,PEP_SD_array);
save('PEP_SD_daily.mat','PEP_SD_array');

LVET_M_array = table2array(LVET_M)';
LVET_M_array = horzcat(subjectCol,LVET_M_array);
save('LVET_M_daily.mat','LVET_M_array');

LVET_SD_array = table2array(LVET_SD)';
LVET_SD_array = horzcat(subjectCol,LVET_SD_array);
save('LVET_SD_daily.mat','LVET_SD_array');

CO_M_array = table2array(CO_M)';
CO_M_array = horzcat(subjectCol,CO_M_array);
save('CO_M_daily.mat','CO_M_array');

CO_SD_array = table2array(CO_SD)';
CO_SD_array = horzcat(subjectCol,CO_SD_array);
save('CO_SD_daily.mat','CO_SD_array');

RMSSD_M_array = table2array(RMSSD_M)';
RMSSD_M_array = horzcat(subjectCol,RMSSD_M_array);
save('RMSSD_M_daily.mat','RMSSD_M_array');

RMSSD_SD_array = table2array(RMSSD_SD)';
RMSSD_SD_array = horzcat(subjectCol,RMSSD_SD_array);
save('RMSSD_SD_daily.mat','RMSSD_SD_array');

RR_M_array = table2array(RR_M)';
RR_M_array = horzcat(subjectCol,RR_M_array);
save('RR_M_daily.mat','RR_M_array');

RR_SD_array = table2array(RR_SD)';
RR_SD_array = horzcat(subjectCol,RR_SD_array);
save('RR_SD_daily.mat','RR_SD_array');