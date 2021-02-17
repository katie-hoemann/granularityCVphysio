%% SVM clustering for ARI data, modified to address reviewers comments
% Cluster with 3 features only: IBI, RSA, PEP
% Notes:
%   - Linear SVM with the function trainSVM_linear.m
%   - All the classifiers across the 10 repetitions are saved in a .mat
%   file at the end
%
% Last modified: Nada Kamona 01/18/2020

%% Prep data:
%  - This is the same as data.mat
clear; clc;

path = cd;
dataPath1 = fullfile(path,'Data','Summary_Physio');
dataPath2 = fullfile(path,'Data','Summary_Sheets');
indexing = 1;
for k =[2 3 6 7 9 11 12 13 14 15 17 18 19 21 22 23 24 27 28 29 32 33 34 36 ...
        37 40 43 45 46 47 49 50 51 54 55 56 57 58 59 61 62 63 64 66 67 68]
    fileName1 =strcat('summary_physio_means_',num2str(k,'%03d'),'.xlsx');
    fileName2 =['Copy of',' ',num2str(k,'%03d'),'_Summary_ES.xlsx'];

    data1 = readtable(fullfile(dataPath1,fileName1));
    All_Subject_Features{indexing}=table2array(data1(:,{'mean_RSA','mean_IBI','mean_PEP'})); 
    All_Subject_Clusters{indexing}=table2array(data1(:,{'C_IDs'}));
    data2 = readtable(fullfile(dataPath2,fileName2));
    AllData{indexing}=string(table2array(data2(:,'EmotionWords'))); 
    indexing =indexing+1;
%AllData{i}= string(table2array(data));
end

%% Clustering
subjectIDs = [2 3 6 7 9 11 12 13 14 15 17 18 19 21 22 23 24 27 28 29 32 33 34 36 ...
        37 40 43 45 46 47 49 50 51 54 55 56 57 58 59 61 62 63 64 66 67 68];
acc = zeros(46,10); % 46 subjects, 1 svm's, 10 reps
accuracy_svm_Linear = zeros(46,1);
trainedClassifiers_all = cell(46, 10); % store all trained classifiers

for sub = 1:46
    % Extract top 3 emotion label words
    A = lower(AllData{1,sub});
    Feat = All_Subject_Features{sub};
    Cluster_ID = All_Subject_Clusters{sub};
    X=[];
    for k = 1:length(A)
        X=[X,strsplit(A(k),',')];
    end
    uniX=unique(X);
    cnt =0;
    count=[];
    for k = uniX
        cnt = cnt +1;
        count(cnt) = sum(X==k);

    end
    maxEmo=string([]);
    for k =1:3
        [val ind]=max(count);
        maxEmo(k)=uniX(ind);
        count(ind)=0;
    end
    Labels = zeros(length(A),1);
    for k = 1:length(A)
       line = strsplit(A(k),',');
       if length(line)>=2
           check = [sum(maxEmo(1)==line), sum(maxEmo(2)==line), sum(maxEmo(3)==line)];
           if sum(check)>=2 % 2+ labels of top emotion words, hence ignore event
               disp(['Skip event #' num2str(k) ' for subject ' num2str(sub)])
               continue;
           end
       end
       
       if sum(maxEmo(1)==line)
          Labels(k)=1;
       elseif sum(maxEmo(2)==line)
           Labels(k)=2;
       elseif sum(maxEmo(3)==line)
           Labels(k)=3;
       end
    end
   
    Feat_sub = Feat(logical(Labels),:);
    Labels_final = Labels(logical(Labels),:);
    Labels_count(sub,:) = [sum(Labels_final==1), sum(Labels_final==2), sum(Labels_final==3)];
    maxEmo_allsubj(sub,:) = maxEmo;
    
    % Prepare for clustering: only use subjects with at least 10 events per
    % emotion label.
    if Labels_count(sub,1)<10 | Labels_count(sub,2)<10 | Labels_count(sub,3)<10
        disp(['Skipping subject ' num2str(sub) ' - Has less than 10 events per emotion label'])
        continue;
    else
        Feat_sub_ordered = [Feat_sub(Labels_final==1,:); Feat_sub(Labels_final==2,:); Feat_sub(Labels_final==3,:)];
        Labels_final_ordered = sort(Labels_final);
        labels_idx{1} = find(Labels_final_ordered==1);
        labels_idx{2} = find(Labels_final_ordered==2);
        labels_idx{3} = find(Labels_final_ordered==3);

        % Training - Repeat 10 times, each with 5 k-fold cross validation
        for rep = 1:10
            disp(['Starting training cycle ' num2str(rep)])
            idx_l1 = randperm(Labels_count(sub,1),10); % label 1, select 10 random events
            idx_l2 = randperm(Labels_count(sub,2),10); % label 2
            idx_l3 = randperm(Labels_count(sub,3),10); % label 3

            % New features set with equal number of events per label
            Feat_sub_4net = [Feat_sub_ordered(labels_idx{1}(idx_l1),:);...
                Feat_sub_ordered(labels_idx{2}(idx_l2),:); ...
                Feat_sub_ordered(labels_idx{3}(idx_l3),:)];
            Labels_final_4net = [ones(10,1); 2*ones(10,1); 3*ones(10,1)];
            
            data_4classifier = [Feat_sub_4net, Labels_final_4net];
            
            % SVM - Linear
            [trainedClassifier, validationAccuracy] = trainSVM_linear(data_4classifier);
            acc(sub, rep) = validationAccuracy;
            trainedClassifiers_all{sub, rep} = trainedClassifier;
        end
    end
    accuracy_svm_Linear(sub,1)=nanmean(acc(sub,:));
end

accuracy_svm_Linear(accuracy_svm_Linear==0,:) = [];
subjectIDs(accuracy_svm_Linear==0,:) = [];

clear line ind indexing k A rep sub path maxEmo uniX X val ...
    validationAccuracy cnt count check Cluster_ID ...
    idx_l1 idx_l2 idx_l3 Feat_sub_4net Labels_final_4net acc Feat ...
    Feat_sub Feat_sub_4net Feat_sub_ordered Labels_final_ordered Labels ...
    Labels_final labels_idx fileName1 fileName2

save('SVMclustering_results_3features_linearSVM.mat')

