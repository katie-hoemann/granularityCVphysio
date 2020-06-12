from __future__ import division
import os
import numpy as np
import math


from utils.ecg_features import preprocess_ecg,ecg_r_peaks
from utils.ecg_features import calculate_instant_ibi,calculate_ibi_rsa
from utils.ecg_quality_check import *

def save_ecg_info(input_file, decimator=1, fs=500):
    s_id = input_file.split('_')[0]
    fs = fs / decimator
    ecg_data = np.loadtxt(input_file, delimiter='\t',usecols=1)

    ecg_data = preprocess_ecg(ecg_data, fs=fs)

    use_samples = int(fs * math.floor(len(ecg_data) / fs))              ##discard samples that don't divide by 500
    seconds = use_samples / fs
    # plot_seconds = int(seconds_per_plot * math.ceil(seconds / seconds_per_plot))                    ##round up to nearest 90 seconds
    ecg_data = ecg_data[:use_samples]
    r_peaks = []
    for i in range(0,len(ecg_data),int(85*fs)):
        r_peaks.extend(i + ecg_r_peaks(ecg_data[i:i+int(85*fs)],fs=fs))
    r_peaks = np.array(r_peaks)
    # ecg_for_30_30 = ecg_data[-120000:-90000]
    # r_peaks = ecg_r_peaks(ecg_for_30_30, fs=fs)
    # r_peaks_temp = r_peaks + 135000
    rpeak_flags = np.zeros(len(ecg_data))
    rpeak_flags[r_peaks] = 1
    # rpeak_flags[r_peaks_temp] = 1
    ibi_b2b = calculate_instant_ibi(r_peaks, fs=fs)
    bad_checks = np.zeros(shape=(len(ecg_data)))
    bad_indices = ibi_quantile_checks(r_peaks, ibi_b2b, neighbors=5,tol=2)
    bad_min_max = ibi_max_min_checks(r_peaks,ibi_b2b,max=2000,min=300)
    if s_id == '061' or s_id == '057' or s_id == '055'\
            or s_id == '050' or s_id == '045' or s_id == '047'\
            or s_id == '036' or s_id == '027' or s_id == '018' or\
            s_id == '013' or s_id == '063' or s_id == '064' or s_id == '012' or \
            s_id == '067' or s_id == '065':
        bad_qc_start,bad_qc_end = ecg_QC_testing(ecg_data,
                                                 r_peaks,
                                                 subchunk_len_flatness=0.8,
                                                 subchunk_len_min_max=1.2,
                                                 subchunk_len_distribution=0.8,
                                                 min_std=1e-5,
                                                 delta_min=-0.005,
                                                 delta_max=0.005,
                                                 feq_sublength=0.064,
                                                 min_skewness=0.25,
                                                 sampling_frequency=fs,
                                                 use_rpeaks=True,
                                                 max_2nd_peak=0.9)
    elif s_id == '061':
        bad_qc_start,bad_qc_end = ecg_QC_testing(ecg_data,
                                                 r_peaks,
                                                 subchunk_len_flatness=0.8,
                                                 subchunk_len_min_max=1.2,
                                                 subchunk_len_distribution=0.8,
                                                 min_std=1e-5,
                                                 delta_min=-0.005,
                                                 delta_max=0.005,
                                                 feq_sublength=0.064,
                                                 min_skewness=0.15,
                                                 sampling_frequency=fs,
                                                 use_rpeaks=True,
                                                 max_2nd_peak=0.99)
    else:
        bad_qc_start,bad_qc_end = ecg_QC_testing(ecg_data,
                                                 r_peaks,
                                                 subchunk_len_flatness=1.2,
                                                 subchunk_len_min_max=1.2,
                                                 subchunk_len_distribution=1.2,
                                                 min_std=1e-5,
                                                 delta_min=-0.005,
                                                 delta_max=0.005,
                                                 feq_sublength=0.064,
                                                 min_skewness=0.45,
                                                 sampling_frequency=fs,
                                                 use_rpeaks=True,
                                                 max_2nd_peak=0.5)

    for locations in bad_indices:
        bad_checks[locations[0]:locations[-1]] = 1

    for locations in bad_min_max:
        bad_checks[locations[0]:locations[-1]] = 2

    for ii in range(len(bad_qc_start)):
        bad_checks[bad_qc_start[ii]:bad_qc_end[ii]] = 3
    ibi = np.interp(np.arange(0, len(ecg_data)), r_peaks[1:], ibi_b2b)
    bad_something = np.zeros_like(bad_checks)
    bad_something[np.where(bad_checks > 0)[0]] = 1
    bad_something[np.where(bad_checks == 0)[0]] = np.nan
    import matplotlib.pyplot as plt
    current_directory = os.getcwd()
    temp_dir = os.path.abspath(os.path.join(current_directory, os.pardir))
    temp_dir = os.path.abspath(os.path.join(temp_dir, os.pardir))
    os.chdir(temp_dir)
    # plots_directory = temp_dir + '/ecg_plots_' + s_id +  '/'
    # if not os.path.isdir(plots_directory):
    #     os.mkdir(plots_directory)
    # # ecg_data = np.zeros(shape=ecg_data.shape)
    # # ecg_data[-120000:-90000] = ecg_for_30_30
    # plt.figure()
    # for plots in range(2):
    #     if plots == 0:
    #         plt.subplot2grid((2,1),(plots,0))
    #         plt.plot(ecg_data[-120000:-90000])
    #         plt.plot(bad_something[-120000:-90000]*ecg_data[-120000:-90000],color='red')
    #         plt.ylabel('ECG')
    #     elif plots == 1:
    #         plt.subplot2grid((2,1),(plots,0))
    #         plt.plot(ibi[-120000:-90000])
    #         plt.plot(bad_something[-120000:-90000]*ibi[-120000:-90000],color='red')
    #         plt.ylabel('IBI')
    # plt.suptitle(input_file.strip('_extended.txt'))
    # plt.savefig(plots_directory + input_file.strip('_extended.txt') + '_ECG.png',dpi=100)
    # plt.close()
    # bad_checks_temp = bad_checks.copy()
    # bad_checks = np.zeros(shape=ecg_data.shape)
    # bad_checks[-120000:-90000] = bad_checks_temp
    # ibi_temp = ibi.copy()
    # ibi = np.zeros(shape=ecg_data.shape)
    # ibi[-120000:-90000] = ibi_temp

    info = np.vstack((ecg_data,bad_checks,rpeak_flags,ibi)).T
    os.chdir(current_directory)
    return info