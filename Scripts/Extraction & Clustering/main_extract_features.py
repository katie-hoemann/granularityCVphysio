import glob
import pandas as pd
import os
from collections import OrderedDict
from utils.helpers import *
from utils.respiration_rate_extraction import calculate_resp_signal

def main_extraction(data_directory,participant_id = None,ECG=True,EDA=False,IMP=False,process_files='*_.txt'):
    if IMP:
        L_info_file = 'Impedance Placement Tracker.xlsx'
        L_info = pd.read_excel(L_info_file)
    participant_directory = data_directory + '/'
    days_directory_list = get_immediate_subdirectories(participant_directory)
    from utils.ecg_parent_functions import save_ecg_info
    from utils.imp_parent_functions import save_imp_info
    for days in days_directory_list:
        if days[:3] != 'Day':
            print ('not a data folder, skipping')
        elif (participant_id == u'002' and days[:5] == 'Day6_') or \
                (participant_id[0:3] == '003' and days[0:5] == 'Day5_') or\
                (participant_id[0:3] == '015' and days[0:5] == 'Day5_') or \
                (participant_id[0:3] == '013' and days[0:5] == 'Day15') or \
                (participant_id[0:3] == '055' and days[0:5] == 'Day2_') or \
                (participant_id[0:3] == '045' and days[0:5] == 'Day14') or \
                (participant_id[0:3] == '050' and days[0:5] == 'Day3_') or \
                (participant_id[0:3] == '056' and days[0:5] == 'Day5_') or \
                (participant_id[0:3] == '014' and days[0:5] == 'Day12') or \
                (participant_id[0:3] == '033' and days[0:5] == 'Day13') or \
                (participant_id[0:3] == '061' and days[0:5] == 'Day13') or \
                (participant_id[0:3] == '062' and days[0:5] == 'Day11') or \
                (participant_id[0:3] == '037' and days[0:5] == 'Day4_') :
            print ('skipping this day for this participant')
        else:
            day_directory = participant_directory + days
            try:
                event_directory = day_directory + '\events/'
                os.chdir(event_directory)
            except:
                pass

            if ECG:
                print (days + ': extracting ecg data and running quality check')

                for file in glob.glob(process_files):
                    ##first column preprocessed ecg, second column quality check, third column r_peaks, 4th column inst_ibi
                    print (file)
                    ecg_info = save_ecg_info(file)
                    current_dir = os.getcwd()
                    if len(ecg_info) >= 60*500:
                        if os.path.isdir(current_dir + '/extracted_info/'):
                            os.chdir(current_dir + '/extracted_info/')
                            np.savetxt(file+'_ECG_info.txt',ecg_info)
                        else:
                            os.mkdir(current_dir + '/extracted_info/')
                            os.chdir(current_dir + '/extracted_info/')
                            np.savetxt(file+'_ECG_info.txt',ecg_info)
                    else:
                        print ('event not long enough')
                    os.chdir(current_dir)
            if IMP:
                print (days + ': extracting impedance data and running quality check')
                for file in glob.glob(process_files):
                    print (file)
                    ##impdata,badpep,badlvet,badensemble,badgrad,instpep,instlvet
                    ids = extract_numbers_from_strings(file)
                    s_id = int(ids[0])
                    day = days[:5]
                    L_info_temp = L_info.at[s_id-1,day]
                    # print L_info
                    imp_data = np.loadtxt(file, delimiter='\t', usecols=3)
                    if len(imp_data) >= 60*500:
                        imp_info = save_imp_info(file,L_info=L_info_temp)
                        current_dir = os.getcwd()
                        if os.path.isdir(current_dir + '/extracted_info/'):
                            os.chdir(current_dir + '/extracted_info/')
                            np.savetxt(file+'_IMP_info.txt',imp_info)
                        else:
                            os.mkdir(current_dir + '/extracted_info/')
                            os.chdir(current_dir + '/extracted_info/')
                            np.savetxt(file+'_IMP_info.txt',imp_info)
                    else:
                        print ('event not long enough')
                    os.chdir(current_dir)


def main_respiration(data_directory,participant_id = None,process_files='*_.txt',
                     baseline_file=r'Baseline_Values.xlsx'):
    participant_directory = data_directory + '/'
    rr_info_filename = participant_directory + 'RR_Info_' + str(participant_id) + '.xlsx'
    days_directory_list = get_immediate_subdirectories(participant_directory)
    baseline_pd = pd.read_excel(baseline_file)
    baseline_pd = baseline_pd[baseline_pd['PPID'] == int(participant_id)]
    baseline_pd_1 = baseline_pd[baseline_pd['Session'] == 1]
    baseline_pd_2 = baseline_pd[baseline_pd['Session'] == 2]
    resp_freq_1 = baseline_pd_1['RespFreq'].values
    resp_freq_2 = baseline_pd_2['RespFreq'].values
    resp_mean = np.zeros(shape=[2,])
    resp_std = np.zeros(shape=[2,])
    resp_mean[0] = np.mean(resp_freq_1)
    resp_mean[1] = np.mean(resp_freq_2)
    resp_std[0] = np.std(resp_freq_1)
    resp_std[1] = np.std(resp_freq_2)
    resp_diff_1 = np.max(np.abs(np.ediff1d(resp_freq_1)))
    resp_diff_2 = np.max(np.abs(np.ediff1d(resp_freq_2)))
    resp_diff_max = np.max([resp_diff_1,resp_diff_2])
    event_counter = 0
    for days in days_directory_list:
        if days[:3] != 'Day':
            print ('not a data folder, skipping')
        elif (participant_id == u'002' and days[:5] == 'Day6_') or \
                (participant_id[0:3] == '003' and days[0:5] == 'Day5_') or\
                (participant_id[0:3] == '015' and days[0:5] == 'Day5_') or \
                (participant_id[0:3] == '013' and days[0:5] == 'Day15') or \
                (participant_id[0:3] == '055' and days[0:5] == 'Day2_') or \
                (participant_id[0:3] == '045' and days[0:5] == 'Day14') or \
                (participant_id[0:3] == '050' and days[0:5] == 'Day3_') or \
                (participant_id[0:3] == '056' and days[0:5] == 'Day5_') or \
                (participant_id[0:3] == '014' and days[0:5] == 'Day12') or \
                (participant_id[0:3] == '033' and days[0:5] == 'Day13') or \
                (participant_id[0:3] == '061' and days[0:5] == 'Day13') or \
                (participant_id[0:3] == '062' and days[0:5] == 'Day11') or \
                (participant_id[0:3] == '037' and days[0:5] == 'Day4_') :
            print ('skipping this day for this participant')
        else:
            day_directory = participant_directory + days
            try:
                event_directory = day_directory + '\events/'
                os.chdir(event_directory)
            except:
                pass

            print (days + ': extracting respiration info')
            for file in glob.glob(process_files):
                ##first column preprocessed ecg, second column quality check, third column r_peaks, 4th column inst_ibi
                print (file)
                ecg_info = event_directory + 'extracted_info/' + file + '_ECG_info.txt'
                imp_info = event_directory + 'extracted_info/' + file + '_IMP_info.txt'
                if not os.path.isfile(ecg_info):
                    print ('event not long enough')
                else:
                    event_counter += 1
                    ecg_info = np.loadtxt(ecg_info)
                    imp_info = np.loadtxt(imp_info)
                    bad_pep = imp_info[:, 1]
                    bad_lvet = imp_info[:, 2]
                    bad_ensemble = imp_info[:, 3]
                    bad_grad = imp_info[:, 4]
                    bad_imp = ((bad_pep + bad_lvet + bad_ensemble + bad_grad) > 0).astype(int)
                    bad_ecg = ecg_info[:,1]
                    resp_info_line = calculate_resp_signal(file,fs=500,figures='no',
                                                      resp_mean=resp_mean,resp_std=resp_std,
                                                      resp_diff=resp_diff_max,bad_ecg=bad_ecg,
                                                      bad_imp=bad_imp,)
                    if event_counter == 1:
                        rr_info_pd = pd.DataFrame(resp_info_line, index=[event_counter, ])
                    else:
                        rr_info_pd = rr_info_pd.append(pd.DataFrame(resp_info_line, index=[event_counter, ]))
    writer = pd.ExcelWriter(rr_info_filename)
    rr_info_pd.to_excel(writer, 'Sheet1')
    writer.save()


def main_extracted_info_to_features(data_directory,participant_id = None,process_files='*_extended.txt'):
    from utils.get_clustering_features import get_summary_features, get_windowed_features
    from utils.ecg_features import preprocess_ecg
    participant_directory = data_directory + '/'
    summary_feature_file = participant_directory +'sitting_features_summary_'+str(participant_id) + '.xlsx'
    rr_file = participant_directory + 'RR_Info_' + str(participant_id) + '.xlsx'
    rr_pd = pd.read_excel(rr_file)
    days_directory_list = get_immediate_subdirectories(participant_directory)
    event_counter = 0
    for days in days_directory_list:
        if days[:3] != 'Day':
            print ('not a data folder, skipping')
        elif (participant_id == u'002' and days[:5] == 'Day6_') or \
                (participant_id[0:3] == '003' and days[0:5] == 'Day5_') or \
                (participant_id[0:3] == '015' and days[0:5] == 'Day5_') or \
                (participant_id[0:3] == '013' and days[0:5] == 'Day15') or \
                (participant_id[0:3] == '055' and days[0:5] == 'Day2_') or \
                (participant_id[0:3] == '045' and days[0:5] == 'Day14') or \
                (participant_id[0:3] == '050' and days[0:5] == 'Day3_') or \
                (participant_id[0:3] == '056' and days[0:5] == 'Day5_') or \
                (participant_id[0:3] == '014' and days[0:5] == 'Day12') or \
                (participant_id[0:3] == '033' and days[0:5] == 'Day13') or \
                (participant_id[0:3] == '061' and days[0:5] == 'Day13') or \
                (participant_id[0:3] == '062' and days[0:5] == 'Day11') or \
                (participant_id[0:3] == '037' and days[0:5] == 'Day4_') :
            print ('skipping this day for this participant')
        else:
            day_directory = participant_directory + days
            try:
                event_directory = day_directory + '\events/'
                extracted_info_directory = event_directory + 'extracted_info/'
                os.chdir(extracted_info_directory)
                os.chdir(event_directory)
            except:
                print ('extracted info folder not found')
                continue
            for file in glob.glob(process_files):
                print (file)
                os.chdir(extracted_info_directory)
                if not os.path.isfile(file + '_ECG_info.txt'):
                    print ("excluding event, less than 60 seconds")
                    continue
                else:
                    ecg_info = np.loadtxt(file + '_ECG_info.txt')
                    imp_info = np.loadtxt(file + '_IMP_info.txt')
                file_rr_pd = rr_pd[rr_pd['Event'] == file]
                event_counter += 1
                sid = int(file.split('_')[0])
                day_id = int(file.split('_')[2])
                event_id = int(file.split('_')[4])
                start_trigger = file.split('_')[6] + '_' +  file.split('_')[7]

                baseline_scl = 0
                baseline_ecg = [0, 0]
                baseline_ibi = baseline_ecg[0]
                baseline_rsa = baseline_ecg[1]

                if len(ecg_info) < 500*60:
                    print ("excluding signal, less than 60 seconds")
                    continue
                else:
                    ecg_info = ecg_info[500*30:]
                    imp_info = imp_info[500*30:]
                ecg_info_original = ecg_info.copy()



                bad_flag,rsa_mean,ibi_mean,ibi_median,pep_mean,lvet_mean,sv_mean,co_mean,rmssd_mean,\
                rsa_std,ibi_std,pep_std,lvet_std,sv_std,co_std,rmssd_std = get_summary_features(ecg_info_original,
                                                                                         ecg_info,
                                                                                         imp_info,
                                                                                         baseline_ibi,
                                                                                         baseline_rsa,)
                bad_ecg = file_rr_pd['Bad_ECG (%)'].values
                bad_imp =  file_rr_pd['Bad_IMP (%)'].values
                bad_any = file_rr_pd['Bad_Any (%)'].values
                bad_all = file_rr_pd['Bad_All (%)'].values
                bad_rr_dif = file_rr_pd['Bad_RR_Dif (%)'].values
                bad_rr_u = file_rr_pd['Bad_RR_U (%)'].values
                bad_rr_l = file_rr_pd['Bad_RR_L (%)'].values
                rr_mean = file_rr_pd['RR_Mean_Clean'].values
                rr_std = file_rr_pd['RR_Std_Clean'].values
                if bad_any > 50.0 or bad_flag:
                    bad_flag = 1
                else:
                    bad_flag = 0
                data_frame_dictionary = OrderedDict([('SID', sid), ('Day', day_id),
                                                     ('Event', event_id), ('Starting Time', start_trigger),
                                                     ('Bad ECG(%)', bad_ecg), ('Bad IMP(%)', bad_imp),
                                                     ('Bad_RR_Dif (%)',bad_rr_dif),('Bad_RR_U (%)',bad_rr_u),
                                                     ('Bad_RR_L (%)',bad_rr_l),
                                                     ('Bad Any(%)', bad_any), ('Bad All(%)', bad_all),
                                                     ('Bad Flag', bad_flag),('Event File',file), ('RSA Mean',rsa_mean),
                                                     ('RSA Std',rsa_std), ('IBI Mean',ibi_mean), ('IBI Std', ibi_std),
                                                     ('PEP Mean',pep_mean), ('PEP Std',pep_std), ('LVET Mean',lvet_mean),
                                                     ('LVET Std',lvet_std), ('CO Mean', co_mean), ('CO Std', co_std),
                                                     ('RMSSD Mean',rmssd_mean), ('RMSSD Std',rmssd_std),
                                                     ('RR Mean',rr_mean),
                                                     ('RR Std',rr_std)])

                if event_counter == 1:
                    data_frame = pd.DataFrame(data_frame_dictionary, index=[event_counter, ])
                else:
                    data_frame = data_frame.append(pd.DataFrame(data_frame_dictionary, index=[event_counter, ]))
    os.chdir(participant_directory)
    writer = pd.ExcelWriter(summary_feature_file)
    data_frame.to_excel(writer, 'Sheet1')
    writer.save()

