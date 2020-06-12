import os
import glob
import numpy as np
import time
import tkinter.filedialog
import pandas as pd
from utils.extraction_utils import load_event

def get_immediate_subdirectories(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]

### start of main code ###

data_to_use = pd.read_excel('IMU_issues_0s_only.xlsx',dtype='object')

drop_data = r"F:\ARI_Data\Sitting_Segments/"  # update to desired destination directory
if not os.path.isdir(drop_data):
    os.makedirs(drop_data)
print('Please Browse to the Data Folder')
time.sleep(2)
data_directory = tkinter.filedialog.askdirectory(initialdir =
                                           r'F:/' )  # server path for data
subject_list = ["{0:03}".format(i) for i in range(19,20)]
count = 0
for participants in subject_list:
    file_directory = data_directory + '/' + participants + '/Phase 2/'
    if os.path.isdir(file_directory):
        print(participants)
        directory_list = get_immediate_subdirectories(file_directory)
        sitting_lengths = []
        days_input = []
        for days in directory_list:
            # print(days)
            if participants in data_to_use.values:
                subject_days = data_to_use.where(data_to_use['SubjectId'] == participants)
                if days[0] == 'D':
                    days_temp = days.replace('y','_')
                    current_day = days_temp.split('_')[1]
                    if current_day in subject_days.values:
                        path = file_directory + days + '/'
                        os.chdir(path)
                        for file in glob.glob("*_event.txt"):
                            if os.stat(file).st_size > 1000:
                                temp_lengths = load_event(event_file_name=file,
                                                          day=days,
                                                          subject_id=participants,
                                                          drop_directory=drop_data,
                                                          threshold=30,count=len(sitting_lengths))
                                sitting_lengths.extend(temp_lengths)
                                days_input.extend(int(current_day)*np.ones(len(temp_lengths)))
                            else:
                                print('ignoring empty file')
        sitting_lengths = np.round(np.array(sitting_lengths),3)
        if count == 0:
            overall_df = pd.DataFrame({participants:sitting_lengths})
            overall_df = pd.concat([overall_df,pd.DataFrame({participants+'_day':days_input})],axis=1)
            count += 1
        else:
            overall_df = pd.concat([overall_df,pd.DataFrame({participants:sitting_lengths})],axis=1)
            overall_df = pd.concat([overall_df,pd.DataFrame({participants+'_day':days_input})],axis=1)

overall_df.to_excel(drop_data + 'sitting_lengths.xlsx')