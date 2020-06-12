######  The followings are helper functions
from __future__ import division
import numpy as np
from scipy.stats import zscore
import math

def ibi_max_min_checks(r_peaks,ibi,max=2000,min=300):
    bad_index = []
    for ii in range(len(ibi)):
        if ibi[ii] >= max or ibi[ii] <= min:
            if ii > 0:
                if ii+2 < len(r_peaks):
                    bad_index.append(r_peaks[ii-1:ii+2])
                else:
                    bad_index.append(r_peaks[ii-1:ii+1])
            else:
                bad_index.append(r_peaks[ii:ii+2])

    return bad_index

def ibi_quantile_checks(r_peaks, ibi,fs=500,neighbors = 10, tol=1):

    ibi_diff = np.abs(np.ediff1d(ibi))
    ibi_bad = np.zeros(shape=len(ibi))
    bad_index = []
    # neighbors = len(ibi_diff)
    if len(ibi_diff)<neighbors:
        neighbors = len(ibi)
    for ii in range(len(ibi_diff)):
        if ii<int(neighbors/2)+1:
            select = np.concatenate((ibi_diff[:ii],ibi_diff[ii+1:neighbors+1]))
            select_ibi = np.concatenate((ibi[:ii],ibi[ii+1:neighbors+1]))
        elif len(ibi_diff)-ii<int(neighbors/2)+1 and len(ibi_diff) - ii > 1:
            select = np.concatenate((ibi_diff[-neighbors-1:ii],ibi_diff[ii+1:]))
            select_ibi = np.concatenate((ibi[-neighbors-1:ii],ibi[ii+1:]))
        elif len(ibi_diff) - ii == 1:
            select = ibi_diff[-neighbors-1:-1]
            select_ibi = ibi[-neighbors-1:-1]
        else:
            select = np.concatenate((ibi_diff[ii - int(neighbors/2):ii],ibi_diff[ii+1:ii+1+int(neighbors/2)]))
            select_ibi = np.concatenate((ibi[ii - int(neighbors/2):ii],ibi[ii+1:ii+1+int(neighbors/2)]))
        qd = quartile_deviation(select)
        med = 3.32*qd
        mad = (np.median(select_ibi) - 2.9*qd)/3
        criterion_beat_diff = (med+mad)/2
        # print 'criterion beat diff: ' + str(criterion_beat_diff)
        if (ibi_diff[ii]) > tol*criterion_beat_diff:
            # print 'ibi difference' + str(ibi_diff[ii])
            if ii+3 < len(r_peaks):
                bad_index.append(r_peaks[ii:ii+4])
            else:
                bad_index.append(r_peaks[ii:ii+3])
            ibi_bad[ii+1] = 1

    return bad_index

def quartile_deviation(data_in):
    q75, q25 = np.percentile(data_in, [75, 25])
    iqr = q75 - q25
    return iqr/2

def ecg_QC_testing(ecg_raw, r_peaks, subchunk_len_flatness=0.8, subchunk_len_min_max=1.2,
                   subchunk_len_distribution=1.2, min_std=1e-5, delta_min=-.005,
                   delta_max=0.005, feq_sublength=0.064,
                   min_skewness=0.45,sampling_frequency=500,
                   use_rpeaks=True,max_2nd_peak=0.5):

    subchunk_len_flatness = subchunk_len_flatness * sampling_frequency
    subchunk_len_min_max = subchunk_len_min_max * sampling_frequency
    subchunk_len_distribution = subchunk_len_distribution * sampling_frequency
    feq_sublength = feq_sublength * sampling_frequency

    bad_start_index=[]
    bad_end_index=[]

############################  Flatness Test
    mm1 = len(ecg_raw)
    kk1 = int(math.floor(mm1 / subchunk_len_flatness))

    for i in range(0, kk1):
        beg = (i) * subchunk_len_flatness + 1
        beg = int(beg)
        ending = (i + 1) * subchunk_len_flatness
        ending = int(ending)
        signal_section = ecg_raw[beg:ending]
        if flat_test(signal_section,min_std)==0:        #Check
            bad_start_index.append(beg)
            bad_end_index.append(ending)

############################  Min_Max Text
    if not use_rpeaks:
        mm2 = len(ecg_raw)
        kk2 = int(math.floor(mm2 / subchunk_len_min_max))

        for i in range(0, kk2):
            beg = (i) * subchunk_len_min_max + 1
            beg = int(beg)
            ending = (i + 1) * subchunk_len_min_max
            ending = int(ending)
            signal_section = ecg_raw[beg:ending]
            if min_max_test(signal_section,delta_min,delta_max) == 0:       #Check
                bad_start_index.append(beg)
                bad_end_index.append(ending)

    else:
        for i in range(0, len(r_peaks) - 1):
            beg = int(r_peaks[i])
            ending = int(r_peaks[i+1]-10)
            signal_section = ecg_raw[beg:ending]
            if len(signal_section) == 0:
                end  = int(r_peaks[i+1])
                signal_section = ecg_raw[beg:end]
            if min_max_test(signal_section,delta_min,delta_max,max_2nd_peak=max_2nd_peak) == 0:       #Check
                bad_start_index.append(r_peaks[i])
                bad_end_index.append(ending)


############################  distribution Text
    if not use_rpeaks:
        mm3 = len(ecg_raw)
        kk3 = int(math.floor(mm3 / subchunk_len_distribution))

        for i in range(0, kk3):
            beg = (i) * subchunk_len_distribution + 1
            beg = int(beg)
            ending = (i + 1) * subchunk_len_distribution
            ending = int(ending)
            signal_section = ecg_raw[beg:ending]
            if len(signal_section) == 0:
                end  = int(r_peaks[i+1])
                signal_section = ecg_raw[beg:end]
            if frequency_test(signal_section,feq_sublength,min_skewness) == 0:        #Check
                bad_start_index.append(beg)
                bad_end_index.append(ending)

########################### use rpeaks

    else:
        for i in range(len(r_peaks)-1):
            if i!=0:
                beg = int(r_peaks[i]-100)
            else:
                beg = int(r_peaks[i])
            ending = int(r_peaks[i+1])-10
            signal_section = ecg_raw[beg:ending]
            if frequency_test(signal_section,feq_sublength,min_skewness) == 0:
                bad_start_index.append(beg)
                bad_end_index.append(ending)


    return (bad_start_index, bad_end_index)

def flat_test(signal_section,min_std):

    check=1
    if np.std(signal_section) < min_std:
        check=0

    return (check)

def min_max_test(signal_section,delta_min,delta_max,max_2nd_peak=0.5):
    import numpy as np

    check=1
    if np.max(signal_section) - np.min(signal_section) < np.max(signal_section):
        check = 0

    if np.mean(signal_section) > 0.5 * np.max(signal_section):
        check = 0

    if np.max(signal_section) - np.min(signal_section) > delta_max:
        check = 0

    if np.max(signal_section) - np.min(signal_section) < delta_min:
        check = 0

    if len(signal_section) > 20:
        if np.max(signal_section[10:-10]) > max_2nd_peak*np.max(signal_section):
            check = 0
            # print ('min_max5')

    return (check)

def frequency_test(signal_section,sublength,min_skewness):
    import math
    import numpy as np
    check = 1
    mm = len(signal_section)

    sum=0
    for i in range(mm):
        sum=sum+abs(signal_section[i])
    abs_avg=sum

    n = int(math.floor(mm / sublength))
    small_sig = []
    for i in range(0,n):

        aaa = signal_section[int((i) * sublength):int((i + 1) * sublength+1)]
        delta = np.max(aaa) - np.min(aaa)
        small_sig.append(delta)

    miu = np.mean(small_sig)
    tau = np.median(small_sig)
    skewness = (miu - tau) / np.mean(abs(small_sig - tau))
    if skewness < min_skewness :
        check = 0
        # print('skewness issue : '+ str(skewness))

    if abs_avg> 0.43:
        check = 0
    return (check)






