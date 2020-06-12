clear;
clc;

%% set parameters
minDuration = 60; % minimum duration (seconds)
generatePlots = 0; % set to 1 to generate plots

%% load & clean data
filename = 'sitting_lengths_0s_only.xlsx';
rawData = importdata(filename);
subjectIDs = rawData.textdata;
data = rawData.data(:,2:end);
missingData = isnan(data(1,:));
data = data(:,~missingData);
subjectIDs = subjectIDs(~missingData);
subjectIDlist = str2double(subjectIDs);

%% get descriptive stats for each participant & generate histogram
for i_subject = 1:length(subjectIDs)
    subjectData = data(:,i_subject);
    subjectData = subjectData(~isnan(subjectData));
    subjectData = subjectData(subjectData>minDuration);
    numRests(i_subject) = numel(subjectData);
    meanRest(i_subject) = mean(subjectData); % duration in seconds
    sdRest(i_subject) = std(subjectData); % duration in seconds
    minRest(i_subject) = min(subjectData); % duration in seconds
    maxRest(i_subject) = max(subjectData); % duration in seconds
    rangeRest(i_subject) = maxRest(i_subject)-minRest(i_subject); % duration in seconds
    totalRest(i_subject) = sum(subjectData); % duration in seconds
    
    subjectData_min = round((subjectData/60),1); % duration in mins
    subjectData_30s_bin = (subjectData_min - mod(subjectData_min,.5))*2; % rests in 30 sec bins
    numBins(i_subject,1) = sum(subjectData_30s_bin); % number of 30 sec bins
    subjectData_1m_bin = subjectData_min - mod(subjectData_min,1); % rests in 1 min bins
    numBins(i_subject,2) = sum(subjectData_1m_bin); % number of 1 min bins
   
%     if generatePlots == 1
    %     figure;
    %     bar1 = bar(numBins(i_subject,:));
    %     barNames = {'30 sec','1 min'};
    %     set(gca,'xticklabel',barNames);
    %     %set(gca,'YScale','log');
    %     ylim([0 4500]);
    %     title(['Number of bins: subject ' char(subjectIDs(i_subject))]);
    %     figureName = [char(subjectIDs(i_subject)) '_num_bins'];
    %     saveas(bar1,figureName,'tiff');
    %     
    %     figure;
    %     edges_min = [.5 1 5 10 30 90]; 
    %     [bins,edges] = discretize(subjectData_min,edges_min);
    %     [counts(i_subject,:),~] = histcounts(subjectData_min,edges_min); % counts for data in mins
    %     catCounts = categorical(bins,[1 2 3 4 5],{'30s-1m','1-5m','5-10m','10-30m','>30m'});
    %     histogram_min = histogram(catCounts)
    %     set(gca,'YScale','log');
    %     ylim([0.95 2500]);
    %     title(['Seated rest duration: subject ' char(subjectIDs(i_subject))]);
    %     figureName = [char(subjectIDs(i_subject)) '_rest_duration'];
    %     saveas(histogram_min,figureName,'tiff');
    %      
    %     close all;
    %end
end

%% convert durations to minutes & generate group-level histograms
meanRest_min = round((meanRest/60),1);
sdRest_min = round((sdRest/60),1);
minRest_min = round((minRest/60),1);
maxRest_min = round((maxRest/60),1);
rangeRest_min = round((rangeRest/60),1);
totalRest_min = round((totalRest/60),1);
%meanCounts = mean(counts,1);

if generatePlots == 1
    figure;
    histogram1 = histogram(meanRest_min);
    ylim([0 22]);
    title('Mean seated rest duration (min)');
    saveas(histogram1,'mean_rest_min','tiff');

    figure;
    histogram2 = histogram(sdRest_min);
    ylim([0 18]);
    title('SD seated rest duration (min)');
    saveas(histogram2,'sd_rest_min','tiff');

    figure;
    histogram3 = histogram(maxRest_min);
    ylim([0 25]);
    title('Max seated rest duration (min)');
    saveas(histogram3,'max_rest_min','tiff');

    figure;
    histogram4 = histogram(numRests);
    ylim([0 18]);
    title('Number of seated rest periods');
    saveas(histogram4,'num_rest','tiff');

    figure;
    histogram5 = histogram(totalRest_min);
    ylim([0 22]);
    title('Total duration of seated rest (min)');
    saveas(histogram5,'total_rest_min','tiff');

    figure;
    bar2 = bar(meanCounts);
    barNames = {'30s-1m','1-5m','5-10m','10-30m','>30m'};
    set(gca,'xticklabel',barNames);
    title('Mean count per duration bin');
    saveas(bar2,'mean_bin_count','tiff');

    figure;
    histogram6 = histogram(counts(:,1));
    ylim([0 18]);
    title('Number of seated rest periods 30s-1m');
    saveas(histogram6,'num_rest_30s-1m','tiff');

    figure;
    histogram7 = histogram(counts(:,2));
    ylim([0 18]);
    title('Number of seated rest periods 1-5m');
    saveas(histogram7,'num_rest_1-5m','tiff');

    figure;
    histogram8 = histogram(counts(:,3));
    ylim([0 18]);
    title('Number of seated rest periods 5-10m');
    saveas(histogram8,'num_rest_5-10m','tiff');

    figure;
    edges = [0 10 20 30 40 50];
    histogram9 = histogram(counts(:,4),'BinEdges',edges);
    ylim([0 20]);
    title('Number of seated rest periods 10-30m');
    saveas(histogram9,'num_rest_10-30m','tiff');

    figure;
    histogram10 = histogram(counts(:,5));
    ylim([0 20]);
    title('Number of seated rest periods >30m');
    saveas(histogram10,'num_rest_30m','tiff');

    figure;
    histogram11 = histogram(numBins(:,1));
    ylim([0 25]);
    title('Number of 30s bins');
    saveas(histogram11,'num_bin_30s','tiff');

    figure;
    histogram12 = histogram(numBins(:,2));
    ylim([0 18]);
    title('Number of 1m bins');
    saveas(histogram12,'num_bin_1m','tiff');

    figure;
    edges = (0:200:5000);
    histogram13 = histogram(numBins(:,1),'BinEdges',edges);
    ylim([0 8]);
    vline([mean(numBins(:,1)) mean(numBins(:,1))-(2*std(numBins(:,1))) mean(numBins(:,1))-std(numBins(:,1))],{'k','r','c'},{'M','-2SD','-1SD'});
    title('Number of 30s bins');
    saveas(histogram13,'num_bin_30s_smaller_bins','tiff');

    figure;
    edges = (0:100:2500);
    histogram14 = histogram(totalRest_min,'BinEdges',edges);
    ylim([0 10]);
    vline([mean(totalRest_min) mean(totalRest_min)-(2*std(totalRest_min)) mean(totalRest_min)-std(totalRest_min)],{'k','r','c'},{'M','-2SD','-1SD'});
    title('Total duration of seated rest (min)');
    saveas(histogram14,'total_rest_min_smaller_bins','tiff');

    close all;
end