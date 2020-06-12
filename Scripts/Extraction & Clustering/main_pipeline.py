from main_extract_features import *
from tkinter import filedialog
import warnings
warnings.filterwarnings("ignore")
import os

print('please browse to the subject folder:')

data_directory = filedialog.askdirectory(initialdir =
                                           r'F:\ARI_Data\Sitting_Segments')

step_1 = True
step_2 = True
step_3 = True
current_directory = os.getcwd()
participant_id = data_directory[-3:]
if step_1:
    print ('converting txt files to extraction info...')
    main_extraction(data_directory,participant_id=participant_id,ECG=True,EDA=False,IMP=True,
                    process_files='*.txt')
    print ('extraction complete!!!')
os.chdir(current_directory)
if step_2:
    print ('do RR calculations')
    main_respiration(data_directory,participant_id=participant_id,process_files='*.txt',
                     baseline_file=r'Baseline_Values.xlsx')
if step_3:
    print ('converting extracted info to usable features...')
    main_extracted_info_to_features(data_directory,participant_id=participant_id,
                                    process_files='*.txt')
    print ('conversion complete!!!')
