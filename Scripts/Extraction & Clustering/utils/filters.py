from __future__ import division
from scipy.signal import ellip, lfilter,filtfilt,savgol_filter,butter


def ellip_bandpass(lowcut, highcut,rs,rp,fs, order=2):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = ellip(order,rs,rp, [low, high], btype='band')
    return b, a

def ellip_bandpass_filter(data, lowcut, highcut,rs,rp, fs = 500, order=2):
    b, a = ellip_bandpass(lowcut, highcut, rs,rp, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y

