from __future__ import division
import numpy as np
import os
import math

from utils.filters import ellip_bandpass_filter
from utils.imp_features import ensemble_average_impedance,get_PEP,get_LVET,get_SV,get_CO
from utils.imp_extraction import get_A_point,get_B_point,get_C_point,get_X_point,get_correct_X,get_correct_B
from utils.imp_quality_check import pep_max_min_checks,lvet_max_min_checks

def save_imp_info(input_file,decimator = 1,
                      fs=500,
                      seconds_per_plot=90,
                      show=False,
                      show_individual=True,
                      features=True,
                      n_ensemble=8,
                      L_info = 15):

    """

    :param input_file:
    :param decimator:
    :param fs:
    :param seconds_per_plot:
    :param show:
    :return:
    """
    s_id = input_file.split('_')[0]
    current_directory = os.getcwd()
    os.chdir(current_directory + '/extracted_info/')
    ecg_info = np.loadtxt(input_file + '_ECG_info.txt')
    ibi = ecg_info[:,3]
    os.chdir(current_directory)
    ecg_bad = ecg_info[:,1]
    ecg_bad = ecg_bad > 0
    RPEAKS = ecg_info[:,2]
    r_peaks = np.where(RPEAKS == 1)[0]
    fs = fs / decimator
    imp_data = np.loadtxt(input_file, delimiter='\t',usecols=3)
    z0 = np.loadtxt(input_file,delimiter='\t',usecols=2)
    imp_data = ellip_bandpass_filter(imp_data, 0.5, 40, rs=1, rp=60, fs=fs, order=2
                                     )

    if len(imp_data) != len(ecg_info):
        imp_data = imp_data[:len(ecg_info)]
        z0 = z0[:len(ecg_info)]
    if len(r_peaks) >= n_ensemble*2:
        averaged_data,bad_ensemble = ensemble_average_impedance(imp_data,r_peaks,ecg_bad,fs,n_ensemble)
    else:
        averaged_data = imp_data #[15000:-15000]
        bad_ensemble = np.ones_like(averaged_data)
        bad_pep = bad_ensemble
        bad_grad = bad_ensemble
        bad_lvet = bad_ensemble
        PEP_full = averaged_data
        LVET_full = averaged_data
        SV_full = averaged_data
        CO_full = averaged_data
        info = np.vstack((averaged_data, bad_pep, bad_lvet, bad_ensemble, bad_grad, PEP_full, LVET_full, SV_full, CO_full))
        return info
    averaged_data = averaged_data
    bad_ensemble = bad_ensemble
    ecg_bad = ecg_bad
    RPEAKS = RPEAKS
    ibi = ibi
    z0 = z0
    r_peaks = np.where(RPEAKS == 1)[0]

    print ('calculating HR')
    n_beats = 61
    if n_beats > len(r_peaks):
        n_beats = len(r_peaks) - 1
    C_ibi = np.zeros(n_beats-1).astype(int)
    r = np.zeros(len(r_peaks)).astype(int)
    for i in range(1,n_beats):
        # print i
        cycle = averaged_data[r_peaks[i]-int(.25*fs):r_peaks[i]+int(.5*fs)]
        r[i] = int(.25 * fs)
        if not cycle.size:
            continue
        C_ibi[i-1] = get_C_point(cycle,r[i]) + r_peaks[i] - r[i]

    beat_length = np.mean(np.ediff1d(C_ibi[:]))
    print ('detecting characteristic points')
    C = np.zeros(len(r_peaks)).astype(int)
    B = np.zeros(len(r_peaks)).astype(int)
    series_B = np.zeros(len(r_peaks)).astype(int)
    series_X = np.zeros(len(r_peaks)).astype(int)
    series_C = np.zeros(len(r_peaks)).astype(int)
    A = np.zeros(len(r_peaks)).astype(int)
    X = np.zeros(len(r_peaks)).astype(int)
    r = np.zeros(len(r_peaks)).astype(int)
    flags_C = np.zeros(len(r_peaks)).astype(int)
    for i in range(len(r_peaks)):
        cycle = averaged_data[r_peaks[i]-int(.25*fs):r_peaks[i]+int(.5*fs)]
        r[i] = int(.25 * fs)
        if not cycle.size:
            flags_C[i] = 1
            continue
        elif len(cycle) < 375:
            flags_C[i] = 1
            continue
        else:
            if i>0:
                C[i] = get_C_point(cycle,r[i],last_C=C[i-1]-r_peaks[i-1])
            else:
                C[i] = get_C_point(cycle,r[i])
            if i == 0:
                b2b_length = beat_length
            else:
                b2b_length = beat_length
            A[i] = get_A_point(cycle,r[i], b2b_length)
            B[i] = get_B_point(cycle, b2b_length, A[i],r[i])
            A[i] = A[i] + r_peaks[i]-int(.25*fs)
            series_B[i] = B[i]
            B[i] = B[i] + r_peaks[i]-int(.25*fs)
            series_C[i] = C[i]
            C[i] = C[i] + r_peaks[i]-int(.25*fs)

    delete_ind = np.where(flags_C == 1)[0]
    r_peaks = np.delete(r_peaks,delete_ind)
    C = np.delete(C, delete_ind)
    B = np.delete(B, delete_ind)
    A = np.delete(A, delete_ind)
    X = np.delete(X, delete_ind)
    r = np.delete(r, delete_ind)
    series_B = np.delete(series_B, delete_ind)
    series_X = np.delete(series_X, delete_ind)
    series_C = np.delete(series_C, delete_ind)

    # print ('correcting B point')
    corrected_B = np.zeros(shape=B.shape)
    corrected_X = np.zeros(shape=X.shape)
    series_B_original = series_B.copy()
    if B[0] == 0:
        corrected_B = get_correct_B(series_B[1:],averaged_data)
        corrected_B = np.append(0,corrected_B)
    else:
        corrected_B = get_correct_B(series_B,averaged_data)

    instantaneous_CC = np.ediff1d(C)

    for i in range(len(r_peaks)):
        cycle = averaged_data[r_peaks[i]-int(.25*fs):r_peaks[i]+int(.5*fs)]
        r[i] = int(.25 * fs)
        if not cycle.size:
            continue
        if i == 0:
            inst_length = instantaneous_CC[1]
        elif i == len(r_peaks) - 1:
            inst_length = instantaneous_CC[len(r_peaks) - 2 ]
        else:
            inst_length = instantaneous_CC[i]
        X[i] = get_X_point(cycle,corrected_B[i],series_C[i],fs,length=inst_length)
        series_X[i] = X[i]
        X[i] = X[i]  +  r_peaks[i]-int(.25*fs)

    series_X_original = series_X.copy()

    # print ('correcting X point')
    if X[0] == 0:
        corrected_X = get_correct_X(series_X[1:],corrected_B,beat_length)
        corrected_X = np.append(0,corrected_X)
    else:
        corrected_X = get_correct_X(series_X,corrected_B,beat_length)

    for i in range(len(r_peaks)):
        cycle = averaged_data[r_peaks[i]-int(.25*fs):r_peaks[i]+int(.5*fs)]
        if not cycle.size:
            continue
        else:
            r[i] = int(.25 * fs)
            corrected_X[i] = corrected_X[i] + r_peaks[i]-int(.25*fs)
            corrected_B[i] = corrected_B[i] + r_peaks[i] - int(.25 * fs)

    X_original = X.copy()
    B_original = B.copy()
    B = corrected_B
    X = corrected_X
    imp_data = averaged_data

    use_samples = int(fs * math.floor(len(averaged_data) / fs))
    seconds = use_samples / fs
    plot_seconds = int(seconds_per_plot * math.ceil(seconds / seconds_per_plot))  ##round up to nearest 90 seconds
    num_subplots = int(plot_seconds / seconds_per_plot)  ##each subplot of 90 s
    ymin = np.max((np.min(averaged_data), -1))
    ymax = np.min((np.max(averaged_data), 2))
    # data_for_plotting = {'Data':averaged_data,'A':A,'B':B_original,'C':C,'X':X_original,
    #                      'Num_Subplots':num_subplots,'ymin':ymin,'ymax':ymax,'seconds_per_plot':seconds_per_plot}
    # dump_file = open('F:\pickle_impedance_files/' + input_file + '_IMP_ensemble.pickle', 'wb')
    # pickle.dump(data_for_plotting, dump_file)
    # dump_file.close()

    print ('calculating LVET and PEP')
    PEP = get_PEP(B,r_peaks,fs)
    PEP_full = np.interp(np.arange(0, len(imp_data)), B, PEP)
    pep_grad = np.gradient(PEP_full)
    bad_pep = pep_max_min_checks(PEP_full,max=200,min=30)

    LVET = get_LVET(X, B, fs)
    SV = get_SV(Z0=z0[C],L=L_info,lvet = LVET,amp_C=imp_data[C])
    SV_full = np.interp(np.arange(0, len(imp_data)), C, SV)
    CO_full = get_CO(SV_full,ibi)
    # print ('mean_SV:' + str(np.mean(SV_full)))
    # print ('mean_CO:' + str(np.mean(CO_full)))
    LVET_full = np.interp(np.arange(0, len(imp_data)), X, LVET)

    lvet_grad = np.gradient(LVET_full)
    bad_lvet = lvet_max_min_checks(LVET_full,max=500,min=100)
    bad_grad = np.zeros(shape=len(LVET_full))
    bad_grad[np.where(np.abs(lvet_grad)>30)[0]] = 1
    bad_grad[np.where(np.abs(pep_grad)>20)[0]] = 2

    bad_something = bad_pep + bad_lvet + bad_ensemble + bad_grad
    bad_something[np.where(bad_something > 0)[0]] = 1
    bad_something[np.where(bad_something == 0)[0]] = np.nan


    z0_premean = np.round(np.mean(z0[-90000:-75000]),3)
    z0_postmean = np.round(np.mean(z0[-75000:-60000]),3)
    temp_dir = os.path.abspath(os.path.join(current_directory, os.pardir))
    temp_dir = os.path.abspath(os.path.join(temp_dir, os.pardir))
    os.chdir(temp_dir)
    z0_file = open('z0_INFO.txt','a+')
    np.savetxt(z0_file,np.reshape(np.vstack((z0_premean,z0_postmean)),newshape=(1,2)),delimiter='\t')
    z0_file.close()
    os.chdir(current_directory)
    current_directory = os.getcwd()
    temp_dir = os.path.abspath(os.path.join(current_directory, os.pardir))
    temp_dir = os.path.abspath(os.path.join(temp_dir, os.pardir))
    os.chdir(temp_dir)

    os.chdir(current_directory)
    info = np.vstack((imp_data,bad_pep,bad_lvet,bad_ensemble,bad_grad,PEP_full,LVET_full,SV_full,CO_full))
    return info.T
