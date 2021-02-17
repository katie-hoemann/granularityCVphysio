ARCHIVE
> Scripts and files that are either outdated or unused in final analyses; stored here as back-up

FUNCTIONS
> MATLAB functions that are necessary for running scripts

INPUT
> Data files that are necessary for running scripts

SCRIPTS
> MATLAB and Python scripts used in final analyses
EXTRACTION & CLUSTERING folder contains Python scripts for data extraction, signal processing, and person-specific clustering (see separate README for details)
seated_rest_merge_subject_data_REVISED.m = compiles and aggregates data from person-specific data extraction and clustering sheets (produced by Python scripts)
seated_rest_merge_subject_data_setNandA_FINAL.m = compiles and aggregates data from person-specific data extraction and clustering sheets (produced by Python scripts) - for analyses using standard prior N and alpha hyperparameters
granularity_calculations_FINAL.m = calculates overall and daily granularity and other measures from emotion intensity ratings in end-of-day surveys
RSA_granularity_correlations_REVISED.m = performs group-level (across-subjects) correlations and regressions examining relationship between overall granularity and resting RSA
RSA_granularity_correlations_withinPerson_FINAL.m = performs participant-level (within-subjects) correlations and regressions examining relationship between daily granularity and resting RSA
clusteringAnalysis_granularity_correlations_REVISED.m = performs group-level (across-subjects) correlations examining relationship between overall granularity and number of clusters discovered in seated rest data
clusteringAnalysis_granularity_ANOVAs_FINAL.m = performs group-level (across-subjects) one-way ANOVAs examining whether granularity differs by participants' assigned cluster (using data from group-level clustering performed using Python script)
main_svm_classification_3features.m = performs person-specific classification analyses using linear SVM
classificationAnalysis_merge_subject_data_REVISED.m = aggregates data from granularity analyses and person-specific classification analyses
classificationAnalysis_granularity_correlations_REVISED.m = performs group-level (across-subjects) correlations examining relationship between overall granularity and classifier performance (on event-based physio data)
grouplevel_clustering_visualization_final.m = performs group-level (across-subjects) visualizations of clustering results
main_clustering_grouplevel_wNMI.py = performs group-level (across-subjects) clustering analyses, including post-hoc assessment of normalized mutual information (NMI) of each feature's contribution

