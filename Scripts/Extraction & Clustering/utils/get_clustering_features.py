from __future__ import division
import numpy as np
from utils.ecg_features import calculate_ibi_rsa,calculate_rmssd

def get_mean_scl_rest(eda_info):
    bad_eda = (eda_info[:,1]).astype(int)
    inst_scl = eda_info[:,2]
    flag = 0
    if np.sum(bad_eda) <= 0.5*(len(bad_eda)):
        mean_scl = np.mean(inst_scl[~bad_eda])
    else:
        flag = 1
        mean_scl = 0

    return flag,mean_scl

def get_mean_ibi_rsa_rest(ecg_info):
    bad_ecg = (ecg_info[:,1]).astype(int)
    r_peaks = ecg_info[:,2]
    inst_ibi = ecg_info[:,3]
    ## calculate rsa
    r_peaks_ind = np.where(r_peaks==1)[0]
    _,_,rsa,bad_rsa_flag = calculate_ibi_rsa(ecg_info[:,0],r_peaks_ind,bad_ecg,window_size=30,overlap=29)
    good_rsa_ind = np.where(bad_rsa_flag==0)[0]
    flag = 0
    if np.sum(bad_ecg) <= 0.5*(len(bad_ecg)):
        mean_ibi = np.mean(inst_ibi[~bad_ecg])
        mean_rsa = np.mean(rsa[good_rsa_ind])
    else:
        flag = 1
        mean_ibi = 0
        mean_rsa = 0

    return flag,mean_ibi,mean_rsa



def get_summary_features(ecg_info_original,ecg_info,imp_info,baseline_ibi,baseline_rsa):

    flag = 0
    bad_ecg = ecg_info[:,1]
    r_peaks = ecg_info[:,2]
    inst_ibi = ecg_info[:,3]
    # ecg_info_original = ecg_info_original[90000:-60000] ##since we don't need the whole, just enough to calculate relevant rsa
    ## calculate rsa
    r_peaks_ind = np.where(r_peaks==1)[0]
    rmssd, bad_rmssd = calculate_rmssd(r_peaks_ind,inst_ibi,bad_ecg)
    _,_,rsa,bad_rsa_flag = calculate_ibi_rsa(ecg_info_original[:,0],
                                             np.where(ecg_info_original[:,2]==1)[0],
                                             ecg_info_original[:,1],
                                             window_size=30,
                                             overlap=29)
    # rsa = rsa[150-before:150+after]
    # bad_rsa_flag = bad_rsa_flag[150-before:150+after]
    good_rsa_ind = np.where(bad_rsa_flag == 0)[0]
    good_rmssd_ind = np.where(bad_rmssd == 0)[0]
    rsa = rsa - baseline_rsa
    inst_ibi = inst_ibi - baseline_ibi

    bad_pep = imp_info[:,1]
    bad_lvet = imp_info[:,2]
    bad_ensemble = imp_info[:,3]
    bad_grad = imp_info[:,4]
    bad_imp = bad_pep + bad_lvet + bad_ensemble + bad_grad
    inst_pep = imp_info[:,5]
    inst_lvet = imp_info[:,6]
    inst_sv = imp_info[:,7]
    inst_co = imp_info[:,8]

    bad_anything = bad_ecg + bad_imp
    bad_anything = bad_anything > 0
    if np.sum(bad_anything) <= 0.5*(len(bad_anything)) and good_rsa_ind.size and ecg_info.size > 0 and good_rmssd_ind.size:
        mean_rsa = np.mean(rsa[good_rsa_ind])
        mean_rmssd = np.mean(rmssd[good_rmssd_ind])
        mean_ibi = np.mean(inst_ibi[~bad_anything])
        median_ibi = np.median(inst_ibi[~bad_anything])
        mean_pep = np.mean(inst_pep[~bad_anything])
        mean_lvet = np.mean(inst_lvet[~bad_anything])
        mean_sv = np.mean(inst_sv[~bad_anything])
        mean_co = np.mean(inst_co[~bad_anything])

        std_rsa = np.std(rsa[good_rsa_ind])
        std_rmssd = np.std(rmssd[good_rmssd_ind])
        std_ibi = np.std(inst_ibi[~bad_anything])
        std_pep = np.std(inst_pep[~bad_anything])
        std_lvet = np.std(inst_lvet[~bad_anything])
        std_co = np.mean(inst_co[~bad_anything])
        std_sv = np.std(inst_sv[~bad_anything])
    else:
        flag = 1
        mean_rsa = 0
        std_rsa = 0
        mean_ibi = 0
        median_ibi = 0
        std_ibi = 0
        mean_pep = 0
        std_pep = 0
        mean_lvet = 0
        std_lvet = 0
        mean_co = 0
        std_co = 0
        mean_sv = 0
        std_sv = 0
        mean_rmssd = 0
        std_rmssd = 0

    return flag,mean_rsa,mean_ibi,median_ibi,mean_pep,mean_lvet,\
           mean_sv,mean_co,mean_rmssd,std_rsa,std_ibi,std_pep,std_lvet,std_sv,std_co,std_rmssd


def get_windowed_features():

    return 0

def calculate_bad_percentages(ecg_info, imp_info):


    bad_ecg = (ecg_info[:,1] > 0).astype(int)
    bad_pep = imp_info[:,1]
    bad_lvet = imp_info[:,2]
    bad_ensemble = imp_info[:,3]
    bad_grad = imp_info[:,4]
    bad_imp = ((bad_pep + bad_lvet + bad_ensemble + bad_grad) > 0).astype(int)
    bad_anything = bad_ecg + bad_imp
    bad_anything = (bad_anything > 0).astype(int)
    bad_everything = bad_ecg * bad_lvet * bad_imp
    bad_everything = (bad_everything > 0).astype(int)

    bad_ecg_percentage = np.round(100*np.sum(bad_ecg)/len(bad_ecg),2)
    bad_imp_percentage = np.round(100*np.sum(bad_imp)/len(bad_imp),2)

    bad_any_percentage = np.round(100*np.sum(bad_anything)/len(bad_anything),2)
    bad_all_percentage = np.round(100*np.sum(bad_everything)/len(bad_everything),2)

    return bad_ecg_percentage,bad_imp_percentage,bad_any_percentage,bad_all_percentage