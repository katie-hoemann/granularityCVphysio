import numpy as np
import matplotlib.pyplot as plt
import math

from scipy.signal import detrend,hamming
# from utils.ecg_extraction import preprocess_ecg,ecg_r_peaks
# from utils.ecg_quality_check import ibi_max_min_checks,ibi_quantile_checks

from biosppy import ecg as ECG
from utils.filters import ellip_bandpass_filter

def preprocess_ecg(ecg,fs=500):

    return ellip_bandpass_filter(ecg, 1, 12, rs=0.15, rp=80, fs=fs, order=2)

def ecg_r_peaks(ecg,fs=500):
    try:
        x = ECG.ecg(ecg,fs,show=False)
        r_peaks = x[2]
    except:
        r_peaks = np.array([-1])
    return r_peaks

def calculate_rmssd(r_peaks_ind,ibi,bad_ecg,window=30,overlap=0):
    rmssd = []
    for i in range(500*window,len(ibi)+1,500*window):
        if i == 500*window:
            r = r_peaks_ind[r_peaks_ind <=i]
        else:
            r = r_peaks_ind[(r_peaks_ind > i - 500*window) & (r_peaks_ind <=i)]
        # relevant_ibi = ibi[r]
        relevant_ibi = []
        flag = 0
        for (ind,values) in enumerate(r):
            if ibi[values] >= 2000:
                relevant_ibi.append(ibi[values]/2)
                relevant_ibi.append(ibi[values]/2)
            elif ibi[values] <= 300 and ind < (len(r) - 1):
                relevant_ibi.append(ibi[values]+ibi[r[ind+1]])
                flag = 1
            elif flag == 0:
                relevant_ibi.append(ibi[values])
            elif flag == 1:
                flag = 0
        relevant_ibi = np.array(relevant_ibi)
        ibi_diff = np.ediff1d(relevant_ibi)
        ibi_diff = ibi_diff ** 2
        ibi_diff = ibi_diff.sum()
        ibi_diff = ibi_diff / len(relevant_ibi)
        rmssd.append(np.sqrt(ibi_diff))
    mean_rmssd = np.mean(rmssd)
    if len(rmssd) > 1:
        std_rmssd = np.std(rmssd)
    bad_rmssd = np.zeros(len(rmssd))
    for (i,r) in enumerate(rmssd):
        if r >= 90 or r <=15:
            bad_rmssd[i] = 1
        else:
            bad_rmssd[i] = 0
    return np.array(rmssd),bad_rmssd

def calculate_instant_ibi(r_peaks,fs=500):
    rr_difference_samples = np.ediff1d(r_peaks)
    return (rr_difference_samples / fs) * 1000

def calculate_ibi_rsa(ecg,r_peaks,bad_ecg=None,fs=500,window_size = 30,overlap = 29):
    ibi_mean = []
    ibi_std = []
    rsa = []
    bad_flag = []
    rr_difference_samples = np.ediff1d(r_peaks)
    ibi_milliseconds = (rr_difference_samples / fs) * 1000
    for ii in range(0,len(ecg),int((window_size - overlap)*fs)):
        if len(ecg[ii:ii+window_size*fs]) >= window_size*fs:
            rpeaks_in_window = ecg_r_peaks(ecg[ii:ii+window_size*fs],fs=fs)
            if (rpeaks_in_window[0] == -1) or (len(rpeaks_in_window)<=2) or (np.sum(bad_ecg[ii:ii+window_size*fs]) >= 0.25*len(ecg[ii:ii+window_size*fs])):
                bad_flag.append(1)
            else:
                bad_flag.append(0)
            # rpeaks_in_window = r_peaks[(r_peaks >= ii) & (r_peaks <= ii + window_size*fs)]
            if rpeaks_in_window[0] == -1 or (len(rpeaks_in_window)<=2):
                ibi_mean.append(-1)
                ibi_std.append(-1)
                rsa.append(-1)
            else:
                rr_difference_samples = np.ediff1d(rpeaks_in_window)
                rr_difference_milliseconds = (rr_difference_samples / fs) * 1000
                ibi_mean.append(np.mean(rr_difference_milliseconds))
                ibi_std.append(np.std(rr_difference_milliseconds))
                rsa.append(calculate_rsa(rr_difference_milliseconds))

    return ibi_mean,ibi_std,np.array(rsa),np.array(bad_flag)

def calculate_rsa(ibi,show_plot=False):
    from scipy.interpolate import CubicSpline
    detrend_ibi = detrend(ibi)
    cubic_coeffs = CubicSpline(np.arange(0,len(detrend_ibi)),detrend_ibi)
    ibi_cubic = cubic_coeffs(np.arange(0,len(detrend_ibi)))
    N = len(ibi_cubic)
    h_window = hamming(N)
    windowed_ibi = np.multiply(detrend_ibi,h_window)
    power_spectrum = np.abs(np.fft.fft(windowed_ibi))**2
    t = np.arange(len(ibi))
    freq = np.fft.fftfreq(t.shape[-1])
    if show_plot:
        plt.plot(freq[:N/2],(power_spectrum[:N/2]))
        plt.show()
    rsa = np.log(np.sum(power_spectrum[(0.12<freq) & (freq<0.4)])*2*(freq[1] - freq[0]))

    return rsa

