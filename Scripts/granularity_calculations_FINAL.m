clear all;
clc;

%% load data file, along with word file that includes raw norms
dataFile = '18ARIEOD_daily_filtered.xlsx';
wordFile = 'words18.csv'; 
rawData = importdata(dataFile);
allData = rawData.data.NoLateSurveys;
subjectIDlist = unique(allData(:,1)); % grab subject IDs from first column of data file
words = readtable(wordFile); 
wordList = rawData.colheaders.NoLateSurveys(5:end)';  % grab sampled words from top row of data file

%% set valence and arousal categories for sampled words
for i_word = 1:height(words) % define valence categories
if words.Valence(i_word) > 5 % derived based on the database mean for Warriner et al (2013)
    valCat(i_word) = {'Positive'};
    positive(i_word) = 1;
    valence(i_word) = 1;
else
    valCat(i_word) = {'Negative'};
    positive(i_word) = 0;
    valence(i_word) = 2;
end
end 

for i_word = 1:height(words) % define arousal categories
if words.Arousal(i_word) > 4.6 % derived based on the sample mean for 88 PANAS-X terms in Warriner et al (2013)
    aroCat(i_word) = {'High'};
    high(i_word) = 1;
else
    aroCat(i_word) = {'Low'};
    high(i_word) = 0;
end
end 

words = [words valCat' aroCat']; % append table with category assignments
words.Properties.VariableNames(5:end) = {'ValCat' 'AroCat'}; % label new variables
labels = [positive' high']; % create matrix for logical indexing in ICC commands

%% grab data for each subject and run through calculations
for i_subject = 1:length(subjectIDlist)
    subjectData = [];
    subjectID = subjectIDlist(i_subject);
    index1 = find(allData(:,1)==subjectID);
    subjectData = allData(index1,:);
    % remove missing data
    missingData = any(isnan(subjectData),2); 
    subjectData = subjectData(~missingData,:); 
    subjectData_ratings = subjectData(:,5:end);
    % calculate granularity (ICCs and zICCs)
    rawICC_subject(i_subject,1) = ICC(subjectData_ratings(:,(labels(:,1)==0)),'A-k'); % negative valence ICC
    rawICC_subject(i_subject,2) = ICC(subjectData_ratings(:,(labels(:,1)==1)),'A-k'); % positive valence ICC
    rawICC_subject(i_subject,3) = (rawICC_subject(i_subject,1)+rawICC_subject(i_subject,2))/2; % valence mean ICC
    rawICC_subject(rawICC_subject<0) = 0;
    rtoZ_subject(i_subject,1) = 0.5*log((1+rawICC_subject(i_subject,1))/(1-rawICC_subject(i_subject,1)));
    rtoZ_subject(i_subject,2) = 0.5*log((1+rawICC_subject(i_subject,2))/(1-rawICC_subject(i_subject,2)));
    rtoZ_subject(i_subject,3) = 0.5*log((1+rawICC_subject(i_subject,3))/(1-rawICC_subject(i_subject,3)));
    % calculate mean and SD for positive and negative valence
    mPositive_subject(i_subject) = mean(mean(subjectData_ratings(:,(valence==1))));
    mNegative_subject(i_subject) = mean(mean(subjectData_ratings(:,(valence==2))));
    sdPositive_subject(i_subject) = mean(std(subjectData_ratings(:,(valence==1))));
    sdNegative_subject(i_subject) = mean(std(subjectData_ratings(:,(valence==2)))); 
    dayIDlist = unique(subjectData(:,2));
    for i_day = 1:length(dayIDlist)
        dayData = [];
        dayID = dayIDlist(i_day);
        index2 = find(subjectData(:,2)==dayID);
        dayData = subjectData(index2,5:end);
        % compute ICCs - positive, negative, valence average; combinations of valence x arousal
        rawICC(i_day,1) = ICC(dayData(:,(labels(:,1)==0)),'A-k'); % negative valence ICC
        rawICC(i_day,2) = ICC(dayData(:,(labels(:,1)==1)),'A-k'); % positive valence ICC
        rawICC(i_day,3) = (rawICC(i_day,1)+rawICC(i_day,2))/2; % valence mean ICC
        rawICC(i_day,4) = ICC(dayData(:,(labels(:,1)==0 & labels(:,2)==0)),'A-k'); % negative valence, low arousal ICC
        rawICC(i_day,5) = ICC(dayData(:,(labels(:,1)==1 & labels(:,2)==1)),'A-k'); % positive valence, high arousal ICC
        rawICC(i_day,6) = ICC(dayData(:,(labels(:,1)==0 & labels(:,2)==1)),'A-k'); % negative valence, high arousal ICC
        rawICC(i_day,7) = ICC(dayData(:,(labels(:,1)==1 & labels(:,2)==0)),'A-k'); % positive valence, low arousal ICC
        % set lower bound of ICCs to 0 (no negative values)
        rawICC(rawICC<0) = 0;
        % calculate granularity from ICCs
        gran(i_day,1) = 1-rawICC(i_day,1);
        gran(i_day,2) = 1-rawICC(i_day,2);
        gran(i_day,3) = 1-rawICC(i_day,3);
        gran(i_day,4) = 1-rawICC(i_day,4);
        gran(i_day,5) = 1-rawICC(i_day,5);
        gran(i_day,6) = 1-rawICC(i_day,6);
        gran(i_day,7) = 1-rawICC(i_day,7);
        % Fisher transform ICCs (noted as zICCs)
        rtoZ(i_day,1) = 0.5*log((1+rawICC(i_day,1))/(1-rawICC(i_day,1)));
        rtoZ(i_day,2) = 0.5*log((1+rawICC(i_day,2))/(1-rawICC(i_day,2)));
        rtoZ(i_day,3) = 0.5*log((1+rawICC(i_day,3))/(1-rawICC(i_day,3)));
        rtoZ(i_day,4) = 0.5*log((1+rawICC(i_day,4))/(1-rawICC(i_day,4)));
        rtoZ(i_day,5) = 0.5*log((1+rawICC(i_day,5))/(1-rawICC(i_day,5)));
        rtoZ(i_day,6) = 0.5*log((1+rawICC(i_day,6))/(1-rawICC(i_day,6)));
        rtoZ(i_day,7) = 0.5*log((1+rawICC(i_day,7))/(1-rawICC(i_day,7)));
        % invert zICCs for intuitive directionality
        zInv(i_day,1) = rtoZ(i_day,1)*-1;
        zInv(i_day,2) = rtoZ(i_day,2)*-1;
        zInv(i_day,3) = rtoZ(i_day,3)*-1;
        zInv(i_day,4) = rtoZ(i_day,4)*-1;
        zInv(i_day,5) = rtoZ(i_day,5)*-1;
        zInv(i_day,6) = rtoZ(i_day,6)*-1;
        zInv(i_day,7) = rtoZ(i_day,7)*-1;
        % calculate mean and SD for positive and negative valence
        mPositive(i_day) = mean(mean(dayData(:,(valence==1))));
        mNegative(i_day) = mean(mean(dayData(:,(valence==2))));
        sdPositive(i_day) = mean(std(dayData(:,(valence==1))));
        sdNegative(i_day) = mean(std(dayData(:,(valence==2)))); 
    end
    %% add subject results to summary table
    ppID = ['PP' num2str(subjectID)];
    if i_subject == 1
        day = (1:1:max(allData(:,2)))';
        day = array2table(day,'VariableNames',{'Day'});
                
        gran_Neg = tablecompile_start(day,dayIDlist,gran(:,1),ppID,'Day','ppDay','ppDay');
        gran_Pos = tablecompile_start(day,dayIDlist,gran(:,2),ppID,'Day','ppDay','ppDay');
        gran_M = tablecompile_start(day,dayIDlist,gran(:,3),ppID,'Day','ppDay','ppDay');
        zInv_Neg = tablecompile_start(day,dayIDlist,zInv(:,1),ppID,'Day','ppDay','ppDay');
        zInv_Pos = tablecompile_start(day,dayIDlist,zInv(:,2),ppID,'Day','ppDay','ppDay');
        zInv_M = tablecompile_start(day,dayIDlist,zInv(:,3),ppID,'Day','ppDay','ppDay');
    else
        gran_Neg = tablecompile_iter(gran_Neg,dayIDlist,gran(:,1),ppID,'Day','ppDay','ppDay');
        gran_Pos = tablecompile_iter(gran_Pos,dayIDlist,gran(:,2),ppID,'Day','ppDay','ppDay');
        gran_M = tablecompile_iter(gran_M,dayIDlist,gran(:,3),ppID,'Day','ppDay','ppDay');
        zInv_Neg = tablecompile_iter(zInv_Neg,dayIDlist,zInv(:,1),ppID,'Day','ppDay','ppDay');
        zInv_Pos = tablecompile_iter(zInv_Pos,dayIDlist,zInv(:,2),ppID,'Day','ppDay','ppDay');
        zInv_M = tablecompile_iter(zInv_M,dayIDlist,zInv(:,3),ppID,'Day','ppDay','ppDay');
    end
    clear rawICC gran rtoZ zInv
end

%% write and save data
% write and save subject-level data
subjectMeasures = horzcat(subjectIDlist, rtoZ_subject, mPositive_subject', mNegative_subject', sdPositive_subject', sdNegative_subject');
variableNames = {'PPID' 'zICC_N','zICC_P','zICC_M','mPos','mNeg','sdPos','sdNeg'};
subjectMeasures_Table = array2table(subjectMeasures,'VariableNames',variableNames);
writetable(subjectMeasures_Table,'granularity_affect_measures.xlsx');
save('granularity_affect_measures.mat','subjectMeasures');

% write day-level data
writetable(gran_Neg,'gran_Neg_daily.xlsx');
writetable(gran_Pos,'gran_Pos_daily.xlsx');
writetable(gran_M,'gran_M_daily.xlsx');
writetable(zInv_Neg,'zInv_Neg_daily.xlsx');
writetable(zInv_Pos,'zInv_Pos_daily.xlsx');
writetable(zInv_M,'zInv_M_daily.xlsx');

% save day-level data
subjectCol = [0; subjectIDlist];

gran_Neg_array = table2array(gran_Neg)';
gran_Neg_array = horzcat(subjectCol,gran_Neg_array);
save('gran_Neg_daily.mat','gran_Neg_array');

gran_Pos_array = table2array(gran_Pos)';
gran_Pos_array = horzcat(subjectCol,gran_Pos_array);
save('gran_Pos_daily.mat','gran_Pos_array');

gran_M_array = table2array(gran_M)';
gran_M_array = horzcat(subjectCol,gran_M_array);
save('gran_M_daily.mat','gran_M_array');

zInv_Neg_array = table2array(zInv_Neg)';
zInv_Neg_array = horzcat(subjectCol,zInv_Neg_array);
save('zInv_Neg_daily.mat','zInv_Neg_array');

zInv_Pos_array = table2array(zInv_Pos)';
zInv_Pos_array = horzcat(subjectCol,zInv_Pos_array);
save('zInv_Pos_daily.mat','zInv_Pos_array');

zInv_M_array = table2array(zInv_M)';
zInv_M_array = horzcat(subjectCol,zInv_M_array);
save('zInv_M_daily.mat','zInv_M_array');