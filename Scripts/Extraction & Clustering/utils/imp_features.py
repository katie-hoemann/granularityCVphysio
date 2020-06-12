from __future__ import division
import numpy as np



def get_SV(Z0,L,lvet,amp_C,rho=135):
    SV = rho * ((L/Z0)**2) * lvet/1000 * np.abs(amp_C)
    return SV

def get_CO(SV,ibi):
    hr = 60*ibi/1000
    CO = (SV/1000)*hr
    return CO
def get_PEP(B,r,fs=500):
    return 1000*(B - r)/fs

def get_LVET(X,B,fs=500):
    return 1000*(X - B)/fs

def ensemble_average_impedance(data,r_peaks,ecg_bad,fs=500,n_ensemble=8):
    averaged_data = data.copy()
    bad_ensemble = np.zeros(len(data))
    for i in range(len(r_peaks)):
        # print i
        if r_peaks[i] + int(.5 * fs) <= len(data):
            beat = np.arange(r_peaks[i] - int(.25 * fs), r_peaks[i] + int(.5 * fs))
        else:
            beat = np.arange(r_peaks[i-1] - int(.25 * fs), r_peaks[i-1] + int(.5 * fs))
        ensemble_window = []
        if i < n_ensemble/2:
            for j in range(n_ensemble):
                if not ecg_bad[r_peaks[j]]:
                    if r_peaks[j] - int(.25 * fs) >=0:
                        ensemble_window.append(data[r_peaks[j] - int(.25 * fs):r_peaks[j] + int(.5 * fs)])
                    else:
                        ensemble_window.append(data[r_peaks[j+1] - int(.25 * fs):r_peaks[j+1] + int(.5 * fs)])
            averaged_data[beat] = np.mean(ensemble_window, axis=0)
            if len(ensemble_window) < n_ensemble/2:
                bad_ensemble[beat] = 1
        elif i >= len(r_peaks) - n_ensemble/2:
            for j in range(1,n_ensemble+1):
                if not ecg_bad[r_peaks[-j]]:
                    data_t = data[r_peaks[-j] - int(.25 * fs):r_peaks[-j] + int(.5 * fs)]
                    if len(data_t) != len(beat):
                        # ensemble_window[j,:] = ensemble_window[j-1,:]
                        data_t = data[r_peaks[-j-1] - int(.25 * fs):r_peaks[-j-1] + int(.5 * fs)]
                        ensemble_window.append(data_t)
                    else:
                        # ensemble_window[j,:] = data_t
                        ensemble_window.append(data_t)
            averaged_data[beat] = np.mean(np.array(ensemble_window), axis=0)
            if len(ensemble_window) < n_ensemble/2:
                bad_ensemble[beat] = 1
        else:
            for j in range(int(n_ensemble/2)):
                if not ecg_bad[r_peaks[i-j]]:
                    ensemble_window.append(data[r_peaks[i-j] - int(.25 * fs):r_peaks[i-j] + int(.5 * fs)])
            for j in range(int(n_ensemble/2),n_ensemble):
                t = j - int(n_ensemble/2) + 1
                if not ecg_bad[r_peaks[i+t]]:
                    data_t = data[r_peaks[i+t] - int(.25 * fs):r_peaks[i+t] + int(.5 * fs)]
                    if len(data_t) == len(beat):
                        ensemble_window.append(data_t)
                    elif len(ensemble_window) > 0:
                        ensemble_window.append(ensemble_window[-1])
            averaged_data[beat] = np.mean(np.array(ensemble_window), axis=0)

            if len(ensemble_window) < n_ensemble/2:
                bad_ensemble[beat] = 1
        if np.sum(np.isnan(averaged_data[beat])) > 0:
            averaged_data[beat] = data[beat]
            bad_ensemble[beat] = 1
    return averaged_data,bad_ensemble

