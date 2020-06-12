from main_extract_features import *
from tkinter import filedialog
import warnings
warnings.filterwarnings("ignore")
import os

print('please browse to the sitting segments data directory folder:')

data_directory = filedialog.askdirectory(initialdir =
                                           r'F:\ARI_Data\Sitting_Segments')
subject_list = ["{0:03}".format(i) for i in range(64,69)] # e.g. (10,13) means 010,011,012
current_directory = os.getcwd()
error_count = 0
for pid in subject_list:
    participant_id = pid
    step_1 = True
    step_2 = True
    step_3 = True
    participant_directory = data_directory + '/' + pid
    os.chdir(current_directory)
    if step_1:
        try:
            print ('converting txt files to extraction info...')
            main_extraction(participant_directory,participant_id=participant_id,ECG=True,EDA=False,IMP=True,
                            process_files='*.txt')
            print ('extraction complete!!!')
        except:
            if error_count == 0:
                os.chdir(current_directory)
                with open("pipeline_ERROR.txt", "w") as text_file:
                    text_file.write('ERROR FOR PARTICIPANT ' + pid + ' in STEP 1' '\n')
                error_count += 1
            else:
                os.chdir(current_directory)
                with open("pipeline_ERROR.txt", "a") as text_file:
                    text_file.write('ERROR FOR PARTICIPANT ' + pid + ' in STEP 1' '\n')

    os.chdir(current_directory)
    if step_2:
        try:
            print ('do RR calculations')
            main_respiration(participant_directory,participant_id=participant_id,process_files='*.txt',
                             baseline_file=r'Baseline_Values.xlsx')
        except:
            if error_count == 0:
                os.chdir(current_directory)
                with open("pipeline_ERROR.txt", "w") as text_file:
                    text_file.write('ERROR FOR PARTICIPANT ' + pid + ' in STEP 2' '\n')
                error_count += 1
            else:
                os.chdir(current_directory)
                with open("pipeline_ERROR.txt", "a") as text_file:
                    text_file.write('ERROR FOR PARTICIPANT ' + pid + ' in STEP 2' '\n')
    if step_3:
        try:
            print ('converting extracted info to usable features...')
            main_extracted_info_to_features(participant_directory,participant_id=participant_id,
                                            process_files='*.txt')
            print ('conversion complete!!!')
        except:
            if error_count == 0:
                os.chdir(current_directory)
                with open("pipeline_ERROR.txt", "w") as text_file:
                    text_file.write('ERROR FOR PARTICIPANT ' + pid + ' in STEP 3' '\n')
                error_count += 1
            else:
                os.chdir(current_directory)
                with open("pipeline_ERROR.txt", "a") as text_file:
                    text_file.write('ERROR FOR PARTICIPANT ' + pid + ' in STEP 3' '\n')


