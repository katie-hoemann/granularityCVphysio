clear;
clc;

%% set parameters
minDuration = 60; % minimum duration (seconds)
dataExcluded = 30; % post-posture change data to exclude (seconds) 
generatePlots = 0; % set to 1 to generate plots

%% load & clean data
filename = 'sitting_lengths_0s_only_with_days.xlsx';
rawData = importdata(filename);
subjectIDs = rawData.textdata;
data = rawData.data(:,2:end);
missingData = isnan(data(1,:));
data = data(:,~missingData);
subjectIDs = subjectIDs(~missingData);
subjectIDs = subjectIDs(1:2:end);
        
numRests_day = nan(15,length(subjectIDs));
totalRest_day = nan(15,length(subjectIDs));
numBins_30s_day = nan(15,length(subjectIDs));
numBins_1m_day = nan(15,length(subjectIDs));
subjects_to_remove = [];
daysUsed = zeros(1,length(subjectIDs));

%% get day-level descriptive stats for each participant & generate histogram
for i_subject = 1:length(subjectIDs)
    subjectData = data(:,(i_subject+i_subject-1));
    subjectData = subjectData(~isnan(subjectData));
  
    subjectDays = data(:,(i_subject+i_subject));
    subjectDays = subjectDays(~isnan(subjectDays));
    uniqueDays = unique(subjectDays);
        
    for i_day = 1:length(uniqueDays)
        dayIndex = find(subjectDays == uniqueDays(i_day));
        dayData = subjectData(dayIndex);
        dayData = dayData(dayData>minDuration);
        dayData = dayData-dataExcluded;
        numRests_day(i_day) = numel(dayData);
        totalRest_day(i_day) = sum(dayData);
        dayData_min = round((dayData/60),1);
        dayData_30s_bin = (dayData_min - mod(dayData_min,.5))*2;
        numBins_30s_day(i_day,i_subject) = sum(dayData_30s_bin);
        numBins_30s_plot(i_day) = sum(dayData_30s_bin);
        dayData_1m_bin = dayData_min - mod(dayData_min,1);
        numBins_1m_day(i_day,i_subject) = sum(dayData_1m_bin);
        numBins_1m_plot(i_day) = sum(dayData_1m_bin);
        numRests(i_day,i_subject) = numel(dayData);
        totalRest(i_day,i_subject) = sum(dayData);
        clear dayIndex dayData dayData_min dayData_30s_bin dayData_1m_bin
    end
    
    daysPresent(i_subject) = numel(find(numBins_30s_plot>0));
%     numBins_30s_plot(numBins_30s_plot<60) = 0;
    
    meanBins_30s(i_subject) = mean(numBins_30s_plot);
    maxBins_30s(i_subject) = max(numBins_30s_plot);
    minBins_30s(i_subject) = min(numBins_30s_plot);
    meanBins_1m(i_subject) = mean(numBins_1m_plot);
    maxBins_1m(i_subject) = max(numBins_1m_plot);
    minBins_1m(i_subject) = min(numBins_1m_plot);
    meanRest(i_subject) = mean(totalRest(:,i_subject));
    meanNumRests(i_subject) = mean(numRests(:,i_subject));
    numBins_30s(i_subject) = sum(numBins_30s_plot);
    numBins_1m(i_subject) = sum(numBins_1m_plot);
 
%     if numBins_30s(i_subject) > 1000
       daysUsed(i_subject) = numel(find(numBins_30s_plot>0));
        
    if generatePlots == 1        
        figure;
        bar1 = bar(numBins_30s_plot);
        barNames = uniqueDays;
        set(gca,'xticklabel',barNames);
        %ylim([0 500]); % overall max 755
        hline(60,'k');
        title(['Number of 30s bins per day: subject ' char(subjectIDs(i_subject))]);
        saveas(bar1,[char(subjectIDs(i_subject)) '_num_30s_bins_day'],'tiff');

        figure;
        bar2 = bar(numBins_1m_plot);
        set(gca,'xticklabel',barNames);
        %ylim([0 250]); % overall max 306
        hline(30,'k');
        title(['Number of 1m bins per day: subject ' char(subjectIDs(i_subject))]);
        saveas(bar2,[char(subjectIDs(i_subject)) '_num_1m_bins_day'],'tiff');
        
        figure;
        bar3 = bar(numRests(:,i_subject));
        set(gca,'xticklabel',barNames);
        title(['Number of rest periods per day: subject ' char(subjectIDs(i_subject))]);
        saveas(bar3,[char(subjectIDs(i_subject)) '_num_rests_day'],'tiff');
    
        close all;
    
%     else
%         subjects_to_remove = [subjects_to_remove subjectIDs(i_subject)];
%     end
    end
    clear subjectDays uniqueDays subjectData numBins_30s_plot numBins_1m_plot totalRest numRests
end

%% generate group-level histograms
daysRemoved = daysPresent-daysUsed;
subjectsRemoved = numel(subjects_to_remove);

if generatePlots == 1
    % figure;
    % histogram1 = histogram(meanBins_1m);
    % ylim([0 30]);
    % ylabel('Number of participants');
    % xlabel('Number of bins');
    % title({'Mean number of 1m bins per day',['min length ' num2str(minRest) ', seconds dropped ' num2str(dataExcluded)]});
    % saveas(histogram1,['mean_num_bins_1m_day_' num2str(minRest) ',' num2str(dataExcluded)],'tiff');
    % 
    % figure;
    % histogram2 = histogram(meanBins_30s);
    % ylim([0 25]);
    % ylabel('Number of participants');
    % xlabel('Number of bins');
    % title({'Mean number of 30s bins per day',['min length ' num2str(minRest) ', seconds dropped ' num2str(dataExcluded)]});
    % saveas(histogram2,['mean_num_bins_30s_day_' num2str(minRest) ',' num2str(dataExcluded)],'tiff');
    % 
    % figure;
    % histogram3 = histogram(minBins_1m);
    % ylim([0 40]);
    % ylabel('Number of participants');
    % xlabel('Number of bins');
    % title({'Minimum number of 1m bins per day',['min length ' num2str(minRest) ', seconds dropped ' num2str(dataExcluded)]});
    % saveas(histogram3,['min_num_bins_1m_day_' num2str(minRest) ',' num2str(dataExcluded)],'tiff');
    % 
    % figure;
    % histogram4 = histogram(minBins_30s);
    % ylim([0 30]);
    % ylabel('Number of participants');
    % xlabel('Number of bins');
    % title({'Minimum number of 30s bins per day',['min length ' num2str(minRest) ', seconds dropped ' num2str(dataExcluded)]});
    % saveas(histogram4,['min_num_bins_30s_day_' num2str(minRest) ',' num2str(dataExcluded)],'tiff');

    figure;
    edges = (50:10:350);
    histogram5 = histogram(meanBins_30s,'BinEdges',edges);
    ylim([0 10]);
    vline([mean(meanBins_30s) mean(meanBins_30s)-(2*std(meanBins_30s)) mean(meanBins_30s)-std(meanBins_30s)],{'k','r','c'},{'M','-2SD','-1SD'});
    ylabel('Number of participants');
    xlabel('Number of bins');
    title({'Mean number of 30s bins per day',['min length ' num2str(minDuration) ', seconds dropped ' num2str(dataExcluded)]});
    saveas(histogram5,['mean_num_bins_30s_day_smaller_bins_' num2str(minDuration) ',' num2str(dataExcluded)],'tiff');

    figure;
    edges = (10:5:120);
    histogram6 = histogram(meanNumRests,'BinEdges',edges);
    ylim([0 30]);
    ylabel('Number of participants');
    xlabel('Number of resting periods');
    title({'Mean number of usable resting periods per day',['min length ' num2str(minDuration) ', seconds dropped ' num2str(dataExcluded)]});
    saveas(histogram6,['mean_num_rests_day_' num2str(minDuration) ',' num2str(dataExcluded)],'tiff');

    %close all;
end