import numpy as np
from utils.filters import ellip_bandpass_filter
import matplotlib.pyplot as plt
from scipy import signal
from collections import OrderedDict

def chunks(lst,n):
    for i in range(0,len(lst),n):
        yield lst[i:i+n]

def process_ordered_peaks(peaks,filename=None,resp_mean=None,resp_std=None,resp_diff=None,
                          bad_ecg_flag=[0],bad_imp_flag=[0],plot=False):
    final_peaks = []
    resp_3sd_u = np.max([resp_mean[0] + 3 * resp_std[0],resp_mean[1] + 3 * resp_std[1]])
    resp_3sd_l = np.max([resp_mean[0] - 3 * resp_std[0],resp_mean[1] - 3 * resp_std[1]])
    final_ecg_flag = []
    final_imp_flag = []
    for (i,list_p) in enumerate(peaks):
        for p in list_p:
            if p > 0.17 and p < .4:
                final_peaks.append(p)
                final_ecg_flag.append(bad_ecg_flag[i])
                final_imp_flag.append(bad_imp_flag[i])
                break
    final_flag_plot = np.array(final_ecg_flag).astype(float)
    final_peaks = np.array(final_peaks)
    final_flag_plot[final_flag_plot == 0] = np.nan
    final_flag_plot = final_flag_plot * final_peaks
    rr_bad_flag_u = np.zeros(len(final_peaks))
    rr_bad_flag_l = np.zeros(len(final_peaks))
    rr_bad_flag_d = np.zeros(len(final_peaks))

    rr_bad_flag_u[np.where(final_peaks > resp_3sd_u)] = 1
    rr_bad_flag_l[np.where(final_peaks < resp_3sd_l)] = 1
    rr_bad_flag_d[np.where(np.abs(np.ediff1d(final_peaks) > resp_diff))[0] + 1] = 1

    if plot:
        fig,ax = plt.subplots(1,2)
        # plt.figure()
        axes = ax[0]
        axes.plot(final_peaks,'x-')
        axes.plot(final_flag_plot,color='red',marker='x')
        axes.axhline(y=np.mean(final_peaks),color='black',label='mean')
        axes.axhline(y=resp_3sd_u,color='green',label='baseline_3sd_u')
        axes.axhline(y=resp_3sd_l,xmin=0,color='red',label='baseline_3sd_l')
        axes.set_ylim((0,.5))
        axes.set_title('Resp_Freq')
        axes.legend()
        axes = ax[1]
        axes.plot(np.abs(np.ediff1d(final_peaks)),'o-')
        axes.axhline(y=resp_diff,label='baseline_max',color='black')
        axes.set_title('Resp_Freq_Change')
        axes.legend()
        fig.tight_layout()
        plt.savefig(filename + '.png')

        plt.close()
    return final_peaks,final_ecg_flag,final_imp_flag,rr_bad_flag_u,rr_bad_flag_l,rr_bad_flag_d
def calculate_resp_signal(input_file,fs=500,figures='save',
                          resp_mean=None,resp_std=None,resp_diff=None,bad_ecg=None,bad_imp=None,
                          ):

    factor = 1
    z0 = np.loadtxt(input_file,delimiter='\t',usecols=2)

    sid = int(input_file.split('_')[0])
    day = int(input_file.split('_')[2])

    if len(z0) >= 60*500:
        z0 = z0[30*500:]
        bad_ecg = bad_ecg[30*500:]
        bad_ecg = bad_ecg > 0
        z0_filt = []
        ordered_peaks = []
        bad_ecg_flag = []
        bad_imp_flag = []
        for (i,imp) in enumerate(list(chunks(z0,30*fs))):
            if len(imp) >= 30*fs:
                bad_ecg_chunk = bad_ecg[i*30*fs : (i+1)*30*fs]
                bad_imp_chunk = bad_imp[i*30*fs : (i+1)*30*fs]

                if np.sum(bad_ecg_chunk) >= 0.5*(len(bad_ecg_chunk)):
                    bad_ecg_flag.append(1)
                else:
                    bad_ecg_flag.append(0)

                if np.sum(bad_imp_chunk) >= 0.5*(len(bad_imp_chunk)):
                    bad_imp_flag.append(1)
                else:
                    bad_imp_flag.append(0)
                window = signal.hanning(len(imp))
                z0_taper = imp * window
                z0_filt_temp = ellip_bandpass_filter(z0_taper, 0.12, .4, rs=1, rp=60, fs=fs/factor, order=2)
                z0_filt_temp = signal.detrend(z0_filt_temp)
                z0_filt_temp = z0_filt_temp - z0_filt_temp.mean()
                z0_filt.extend(z0_filt_temp)

                power_spectrum = np.abs(np.fft.fft(z0_filt_temp)) **2

                freq = np.fft.fftfreq(z0_filt_temp.size,(1/(fs/factor)))
                idx = np.where((freq < 0.5) & (freq >= 0))[0]
                peak_frequencies = freq[idx][signal.find_peaks(power_spectrum[idx])[0]]
                temp_ordered_peaks = peak_frequencies[np.argsort(power_spectrum[idx][signal.find_peaks(power_spectrum[idx])[0]])[::-1]]
                ordered_peaks.append(temp_ordered_peaks)
                if figures == 'show':
                    plt.figure()
                    plt.plot(z0_filt_temp)
                    plt.title('dz_filt, chunk ' + str(i))
                    plt.show()
                    plt.figure()
                    plt.plot(freq[idx], (power_spectrum[idx]))
                    plt.title('power_spectrum, chunk ' + str(i))
                    plt.show()
                elif figures == 'save':
                    plt.figure()
                    plt.plot(freq[idx], (power_spectrum[idx]))
                    plt.title('power_spectrum, chunk ' + str(i))
                    plt.savefig()
        final_peaks,final_ecg_flag,final_imp_flag,rr_bad_flag_u,rr_bad_flag_l,rr_bad_flag_d =\
            process_ordered_peaks(ordered_peaks,input_file,resp_mean,resp_std,resp_diff,
                                  bad_ecg_flag=bad_ecg_flag,bad_imp_flag=bad_imp_flag,plot=False)
        bad_any = final_ecg_flag + rr_bad_flag_u + rr_bad_flag_l + rr_bad_flag_d + final_imp_flag
        bad_any = bad_any > 0
        bad_all = final_ecg_flag * rr_bad_flag_u * rr_bad_flag_l * rr_bad_flag_d * final_imp_flag
        bad_all = bad_all > 0
        final_peaks_post = final_peaks[~bad_any]
        if len(final_peaks) > 0:
            data_dictionary = OrderedDict([('SID', sid), ('Day', day),
                                          ('Bad_All (%)', np.round(np.sum(bad_all) / len(final_peaks) * 100, 2)),
                                          ('Bad_Any (%)', np.round(np.sum(bad_any) / len(final_peaks) * 100, 2)),
                                          ('Bad_ECG (%)', np.round(np.sum(final_ecg_flag) / len(final_peaks) * 100, 2)),
                                          ('Bad_IMP (%)', np.round(np.sum(final_imp_flag) / len(final_peaks) * 100, 2)),
                                          ('Bad_RR_Dif (%)', np.round(np.sum(rr_bad_flag_d) / len(final_peaks) * 100, 2)),
                                          ('Bad_RR_L (%)', np.round(np.sum(rr_bad_flag_l) / len(final_peaks) * 100, 2)),
                                          ('Bad_RR_U (%)', np.round(np.sum(rr_bad_flag_u) / len(final_peaks) * 100, 2)),
                                          ('Event', input_file),
                                          ('RR_Mean_Clean', np.mean(final_peaks_post)),
                                           ('RR_Mean_Total', np.mean(final_peaks)),
                              ('RR_Std_Total', np.std(final_peaks)),
                              ('RR_Std_Clean', np.std(final_peaks_post)),
                              ('Bins Left (%)', np.round(len(final_peaks_post) / len(final_peaks) * 100, 2)),
                              ('Total Bins', len(final_peaks))])
        else:
            print ('not peaks in range')
            data_dictionary = OrderedDict([('SID', sid), ('Day', day), ('Bad_All (%)', 100.0),
                                           ('Bad_Any (%)', 100.0),
                                           ('Bad_ECG (%)', 100.0),
                                           ('Bad_RR_Dif (%)', 100.0),
                                           ('Bad_RR_L (%)', 100.0),
                                           ('Bad_RR_U (%)', 100.0),
                                           ('Event', input_file),
                                           ('RR_Mean_Clean', -1),
                                           ('RR_Std_Total', -1),
                                           ('RR_Std_Clean',-1),
                                           ('RR_Mean_Total', -1),
                                           ('Bins Left (%)', -1),
                                           ('Total Bins', 0)])

    else:
        print('event not long enough')
    return data_dictionary
