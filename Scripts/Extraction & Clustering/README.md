# Sitting_Analysis
Please run scripts in the following order:

1. main_extract_sitting_info.py:
This script takes in the raw mindware generated text files for all participants in the range specified in line 24.
It extracts all sitting segments that did not have a prompt and converts saves text files for these in subject folders in the
directory specified in line 17 of the code.

2. main_pipeline.py:
This script runs on a per subject level. When prompted navigate to the subject specific folder created from the last step.
step 1, 2 and 3 should all be kept True for all steps to complete together. If there's an error in later steps we can re-run just 
that step. Once it says "conversion complete", you will have an excel sheet with all the features and quality check info in the 
subject folder with the name "sitting_feautres_summary_[subject id].xlsx

3. main_clustering.py:
This script also runs on a per subject level. When prompted navigate to the subject specific folder containing the excel sheet from last step. Once completed this will create another excel sheet with all the information from the input sheet plus clustering results, the sheet will be named "sitting_features_summary_results_[subject id].xlsx

OPTIONAL: 4. main_clustering_setNandA.py:
This script is identical to the above, except that it implements (user-specific) standard values for the prior N components and alpha hyperparameters.

OPTIONAL: 5. main_marginal_alpha_N.py:
This script runs on a per subject level to perform a Bayesian hyperparameter selection approach known as Empirical Bayes, which examines the  marginal likelihood (i.e., evidence) of model fit across a distribution of possible values for DP-GMM hyperparameters for prior N components and alpha.
