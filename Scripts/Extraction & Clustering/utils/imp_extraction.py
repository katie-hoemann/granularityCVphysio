from __future__ import division
import numpy as np
from scipy.signal import argrelextrema
from scipy.signal import butter, filtfilt
from statsmodels.robust.scale import mad
import biosppy
from spectrum import arburg

from utils.helpers import zero_runs


def get_C_point(beat, R,last_C=0):  # dzdtmax
    input_beat = beat.copy()
    beat = input_beat[R:300]
    local_max_ind = argrelextrema(beat, np.greater, order=1)[0]
    local_max_values = np.argsort(beat[local_max_ind])[::-1]
    bad_flag = 0
    if len(local_max_values) > 1:
        potential_Cs = local_max_ind[local_max_values[:2]]
        C = potential_Cs[0]
        if potential_Cs[0] < potential_Cs[1]:
            C = potential_Cs[0]
            if C < 25:
                C = potential_Cs[1]
        else:
            if potential_Cs[0] - potential_Cs[1] < 20:
                C = potential_Cs[0]
                if C < 25:
                    C = potential_Cs[1]
            elif 100 * (beat[potential_Cs[0]] - beat[potential_Cs[1]]) / beat[potential_Cs[0]] < 5:
                C = potential_Cs[1]
                if C < 25:
                    C = potential_Cs[0]
            else:
                C = potential_Cs[0]
                if C < 25:
                    C = potential_Cs[1]
    elif len(local_max_ind) == 1:
        C = local_max_ind[0]
    else:
        beat = input_beat[R:]
        local_max_ind = argrelextrema(beat, np.greater, order=1)[0]
        local_max_values = np.argsort(beat[local_max_ind])[::-1]
        if len(local_max_values) > 1:
            potential_Cs = local_max_ind[local_max_values[:2]]
            C = potential_Cs[0]
            if potential_Cs[0] < potential_Cs[1]:
                C = potential_Cs[0]
                if C < 25:
                    C = potential_Cs[1]
            else:
                if potential_Cs[0] - potential_Cs[1] < 20:
                    C = potential_Cs[0]
                    if C < 25:
                        C = potential_Cs[1]
                elif 100 * (beat[potential_Cs[0]] - beat[potential_Cs[1]]) / beat[potential_Cs[0]] < 5:
                    C = potential_Cs[1]
                    if C < 25:
                        C = potential_Cs[0]
                else:
                    C = potential_Cs[0]
                    if C < 25:
                        C = potential_Cs[1]
        elif len(local_max_ind) == 1:
            C = local_max_ind[0]
        else:
            bad_flag = 1
            # raise ValueError ('something wrong with C')
    if bad_flag:
        return last_C+R
    else:
        return C + R


def get_A_point(beat, R, length):
    # length = len(beat)
    C = get_C_point(beat, R)
    if length / 3 > C:
        length = len(beat)
    data_look = beat[int(C - length / 3):C]
    if data_look.size:
        data_look = biosppy.tools.smoother(signal=data_look, kernel='bartlett')[0]
    else:
        data_look = beat[:C]
        data_look = biosppy.tools.smoother(signal=data_look, kernel='bartlett')[0]
    A_candidates = argrelextrema(data_look, np.less)[0]
    if A_candidates.size:
        if len(A_candidates) > 1:
            A = A_candidates[np.argmin(data_look[A_candidates])]
            A_star = A_candidates[-1]
            if data_look[A_star] - data_look[A] < 0.1 * data_look[A]:
                A = A_star
        else:
            A = A_candidates
    else:
        A = 0
    return A + int(C - length / 3)


def get_dzdt_H(beat, A, R):
    C = get_C_point(beat, R)
    return beat[C] - beat[A]


def get_d2zdt2(beat, fs=500):
    gradient = np.gradient(beat)
    return gradient


def get_d3zdt3(beat, fs=500):
    dz2dt2 = np.gradient(beat)
    gradient = np.gradient(dz2dt2)
    return gradient


def monotone_segments(A2C):
    d2zdt2 = get_d2zdt2(A2C)
    flag = np.zeros(len(A2C))
    segments = {}
    count = 0
    for i in range(len(A2C)):
        if flag[i] != 1 and d2zdt2[i] >= 0:
            segments[count] = []
            segments[count].append(i)
            flag[i] = 1
            for j in range(i + 1, len(A2C)):
                if d2zdt2[j] >= 0:
                    segments[count].append(j)
                    flag[j] = 1
                else:
                    break
            count = count + 1
    return segments


def get_prominent_segment(beat, A, R):
    C = get_C_point(beat, R)
    A = R
    if len(beat[A:C]) < 3:
        #print ('R and C too close')
        return []
    positive_segments = monotone_segments(beat[A:C])
    prominent_segs = []
    heights = []
    count = 0
    flag = 0
    label = np.zeros(len(positive_segments))
    seg_ind = []
    for i in range(len(positive_segments)):
        seg = positive_segments[i]
        test = beat[seg[0] + A] <= 0.5 * beat[C] and beat[seg[-1] + A] >= (2. / 3) * beat[C]
        if test:
            prominent_segs.append(seg)
            heights.append(beat[A:C][seg[-1]])
            seg_ind.append(i)
            count += 1
    if not heights:
        for i in range(len(positive_segments)):
            seg = positive_segments[i]
            if beat[seg[0] + A] <= 0.5 * beat[C]:
                prominent_segs.append(seg)
                heights.append(beat[A:C][seg[-1]])
                seg_ind.append(i)
                count += 1
    if not heights:
        for i in range(len(positive_segments)):
            seg = positive_segments[i]
            if beat[seg[0] + A] <= 0.75 * beat[C]:
                prominent_segs.append(seg)
                heights.append(beat[A:C][seg[-1]])
                seg_ind.append(i)
                count += 1
        # raise ValueError('No prominent segment')
    if not heights:
        for i in range(len(positive_segments)):
            seg = positive_segments[i]
            if beat[seg[0] + A] <= 0.75 * beat[C]:
                prominent_segs.append(seg)
                heights.append(beat[A:C][seg[-1]])
                seg_ind.append(i)
                count += 1
    if not heights:
        if len(positive_segments) == 0:
            #print ('no positive segments')
            return []
        else:
            return positive_segments[0] + A
    else:
        most_prominent = seg_ind[np.argmax(heights)]
        return positive_segments[most_prominent] + A


def get_inflection_points(beat, A, R, fs=500):
    most_prominent_seg_index = get_prominent_segment(beat, A, R)
    if len(most_prominent_seg_index) == 0:
        analyze_segment = []
    else:
        most_prominent_seg = beat[most_prominent_seg_index]
        H = get_dzdt_H(beat, A, R)
        zero_crossing_kept = []
        analyze_segment = most_prominent_seg[:int(len(most_prominent_seg) / 1.25)]
    if len(analyze_segment) < 3:
        #print ('segment not long enough')
        answer = []
    else:
        d2_seg = get_d2zdt2(analyze_segment, fs)
        d3_seg = get_d3zdt3(analyze_segment, fs)
        zero_crossing_detected = biosppy.tools.zero_cross(signal=d3_seg, detrend=False)[0]
        if not zero_crossing_detected.size:
            #print('no zero crossing detected')
            answer = []
        else:
            for i in zero_crossing_detected:
                if 0 < d2_seg[i] <= 2 * H / fs:
                    zero_crossing_kept.append(i)
            if not zero_crossing_kept:
                for i in zero_crossing_detected:
                    if 0 < d2_seg[i] <= 4 * H / fs:
                        zero_crossing_kept.append(i)
            if not zero_crossing_kept:
                for i in zero_crossing_detected:
                    if d2_seg[i] <= 6 * H / fs:
                        zero_crossing_kept.append(i)
            answer = zero_crossing_kept + most_prominent_seg_index[0]
    return np.array(answer), most_prominent_seg_index


def get_local_max(beat, A, R, most_prominent_seg_index, fs=500):
    most_prominent_seg = beat[most_prominent_seg_index]
    H = get_dzdt_H(beat, A, R)
    analyze_segment = most_prominent_seg[:int(len(most_prominent_seg) / 1.25)]
    if len(analyze_segment) < 3:
        #print ('segment not long enough')
        answer = []
    else:
        d3_seg = get_d3zdt3(analyze_segment, fs)
        local_max = argrelextrema(d3_seg, np.greater)[0]
        if not local_max.size:
            answer = []
        else:
            accept_local_max = np.where(d3_seg[local_max] > 7 * (H / fs))[0]
            answer = local_max[accept_local_max] + most_prominent_seg_index[0]

    return np.array(answer)


def get_inflection_points_relaxed(beat, A, R, B_old, fs=500):
    most_prominent_seg_index = get_prominent_segment(beat, A, R)
    if len(most_prominent_seg_index) == 0:
        analyze_segment = []
    else:
        most_prominent_seg = beat[most_prominent_seg_index]
        B_old = B_old - most_prominent_seg_index[0] + 1
        H = get_dzdt_H(beat, A, R)
        zero_crossing_kept = []
        analyze_segment = most_prominent_seg[B_old:int(len(most_prominent_seg) / 1.25)]
    if len(analyze_segment) < 3:
        #print ('segment not long enough')
        answer = []
    else:
        d2_seg = get_d2zdt2(analyze_segment, fs)
        d3_seg = get_d3zdt3(analyze_segment, fs)
        zero_crossing_detected = biosppy.tools.zero_cross(signal=d3_seg, detrend=False)[0]
        if not zero_crossing_detected.size:
            #print('no zero crossing detected')
            answer = []
        else:
            for i in zero_crossing_detected:
                if d2_seg[i] > 0 and d2_seg[i] <= 2 * H / fs:
                    # if (2 * (H / fs) - d2_seg[i]) / (6 * (H / fs)) > 0.3:
                    zero_crossing_kept.append(i)
            if not zero_crossing_kept:
                for i in zero_crossing_detected:
                    if d2_seg[i] > 0 and d2_seg[i] <= 4 * H / fs:
                        zero_crossing_kept.append(i)
            if not zero_crossing_kept:
                for i in zero_crossing_detected:
                    if d2_seg[i] <= 6 * H / fs:
                        zero_crossing_kept.append(i)

            if not zero_crossing_kept:
                for i in zero_crossing_detected:
                    if d2_seg[i] <= 8 * H / fs:
                        zero_crossing_kept.append(i)

            if not zero_crossing_kept:
                for i in zero_crossing_detected:
                    if d2_seg[i] <= 10 * H / fs:
                        zero_crossing_kept.append(i)

            flag_no_zero = 0
            if not zero_crossing_kept:
                zero_crossing_kept.append(zero_crossing_detected[np.argmin(d2_seg[zero_crossing_detected])])
                flag_no_zero = 1
            answer = zero_crossing_kept + most_prominent_seg_index[0] + B_old
    return np.array(answer), most_prominent_seg_index


def get_local_max_relaxed(beat, A, R, most_prominent_seg_index, B_old, fs=500):
    most_prominent_seg = beat[most_prominent_seg_index]
    if len(most_prominent_seg_index) == 0:
        analyze_segment = []
    else:
        B_old = B_old - most_prominent_seg_index[0] + 1
        H = get_dzdt_H(beat, A, R)
        analyze_segment = most_prominent_seg[B_old:int(len(most_prominent_seg) / 1.25)]
    if len(analyze_segment) < 3:
        #print ('segment not long enough')
        answer = []
    else:
        d3_seg = get_d3zdt3(analyze_segment, fs)
        local_max = argrelextrema(d3_seg, np.greater)[0]
        accept_local_max = []
        if not local_max.size:
            #print('no rapid slope change detected')
            answer = []
        else:
            accept_local_max = np.where(d3_seg[local_max] > 7 * (H / fs))[0]
            flag_no_local_max = 0
            if not accept_local_max:
                accept_local_max = []
                accept_local_max.append(np.argmax(d3_seg[local_max]))
                flag_no_local_max = 1
            answer = local_max[accept_local_max] + B_old + most_prominent_seg_index[0]

    return np.array(answer)

def get_B_point(beat, length, A, R, fs=500, min_RB=int(30 / 2)):
    flag = 0
    zero_crossing, most_prom_seg_ind = get_inflection_points(beat, A, R, fs)
    if not flag:
        # C = get_C_point(beat, R)
        local_max = get_local_max(beat, A, R, most_prom_seg_ind, fs)
        if zero_crossing.size and local_max.size:
            if zero_crossing[-1] >= local_max[-1]:
                B = zero_crossing[-1]
            else:
                B = local_max[-1]
        elif zero_crossing.size and not local_max.size:
            B = zero_crossing[-1]

        elif local_max.size and not zero_crossing.size:
            B = local_max[-1]

        elif not local_max.size and not zero_crossing.size:
            if len(most_prom_seg_ind) > 0:
                B = most_prom_seg_ind[0]
            else:
                B = R

        if B < A:
            B = A
        if B - R < min_RB:
            B_old = R + min_RB
            zero_crossing, most_prom_seg_ind = get_inflection_points_relaxed(beat, A, R, B_old, fs)
            if not flag:
                C = get_C_point(beat, R)
                local_max = get_local_max_relaxed(beat, A, R, most_prom_seg_ind, B_old, fs)
                if zero_crossing.size and local_max.size:
                    if zero_crossing[-1] >= local_max[-1]:
                        B = zero_crossing[-1]
                    else:
                        B = local_max[-1]
                elif zero_crossing.size and not local_max.size:
                    B = zero_crossing[-1]

                elif local_max.size and not zero_crossing.size:
                    B = local_max[-1]

                elif not local_max.size and not zero_crossing.size:
                    B = R + min_RB

                if B < A:
                    B = A
            flag_min = 1

    return B


def find_outliers(B):
    median = np.median(B)
    smad = mad(B)
    outliers_B = np.where(np.abs(B - median) >= 3 * smad)[0]
    return outliers_B

def get_correct_B(B, data, fs=500):
    starting_B = B.copy()

    b, a = butter(4, 0.1)
    B_baseline = filtfilt(b, a, B)
    B_stationary = B - B_baseline

    outliers = find_outliers(B)
    B_stationary[outliers] = np.nan
    B_stationary_no_outliers = np.isnan(B_stationary).astype(int)
    no_outlier_sequences = zero_runs(B_stationary_no_outliers)

    lst = no_outlier_sequences
    temp_list = [x for x in lst if x[1] - x[0] >= 30]
    no_outlier_sequences = temp_list
    n_outliers = np.zeros(shape=(20,))
    n_outliers[1] = 1
    count = 0
    count_iter = 0
    max_iter = 20
    while len(outliers) > 0:
        # #print(len(outliers))
        #print count_iter
        if len(np.unique(n_outliers)) == 1 or count_iter>=max_iter:
            break
        if count < len(n_outliers):
            n_outliers[count] = len(outliers)
            count += 1
        else:
            count = 0
        for o in outliers:
            predicted = 0
            forward_seq = [x for x in no_outlier_sequences if x[1] <= o]
            if forward_seq:
                temp_stationary = B_stationary.copy()
                forward_seq = forward_seq[-1]  # pick the closest sequence to the outlier
                forward_seq_data = B_stationary[forward_seq[0]:forward_seq[1]]
                try:
                    AR_f, p, f = arburg(forward_seq_data, order=len(forward_seq_data) - 1, criteria='AICc')
                except:
                    break
                m = len(AR_f)
                predicted_f = 0.0
                for i in range(forward_seq[1], o + 1):
                    for j in range(m):
                        predicted_f += AR_f[j] * temp_stationary[i - 1 - j]
                    temp_stationary[i] = predicted_f
                predicted_f = predicted_f + B_baseline[o]
            else:
                predicted_f = 0

            backward_seq = [x for x in no_outlier_sequences if x[0] > o]
            if backward_seq:
                temp_stationary = B_stationary.copy()
                backward_seq = backward_seq[0]
                backward_seq_data = B_stationary[backward_seq[0]:backward_seq[1]]
                try:
                    AR_b, p, f = arburg(backward_seq_data, order=len(backward_seq_data) - 1, criteria='AICc')
                except:
                    break
                m = len(AR_b)
                predicted_b = 0.0
                for i in range(backward_seq[0] - 1, o - 1, -1):
                    for j in range(m):
                        predicted_b += AR_b[j] * temp_stationary[i + 1 + j]
                    temp_stationary[i] = predicted_b
                predicted_b = predicted_b + B_baseline[o]
            else:
                predicted_b = 0

            if predicted_b != 0 and predicted_f != 0:
                predicted = (predicted_b + predicted_f) / 2
            elif predicted_b == 0:
                predicted = predicted_f
            elif predicted_f == 0:
                predicted = predicted_b
            # else:
            #     #print ('clean up in aisle outlier!!')
            if np.abs(predicted) >= 200:  # and np.abs(predicted) > starting_B[o]:
                B[o] = starting_B[o]
            elif not np.isnan(predicted):
                B[o] = int(round(predicted))

        B_baseline = filtfilt(b, a, B)
        B_stationary = B - B_baseline
        outliers = find_outliers(B)
        B_stationary[outliers] = np.nan
        B_stationary_no_outliers = np.isnan(B_stationary).astype(int)
        no_outlier_sequences = zero_runs(B_stationary_no_outliers)

        ## remove sequences that are too short
        lst = no_outlier_sequences
        temp_list = [x for x in lst if x[1] - x[0] >= 30]
        no_outlier_sequences = temp_list
        count_iter+=1
    flag = 0
    for i in range(len(B)):
        if B[i] < 140:
            B[i] = starting_B[i]
            flag = 1
        elif np.abs(B[i] - np.mean(starting_B)) > np.abs(starting_B[i] - np.mean(starting_B)):
            B[i] = starting_B[i]
            flag = 1

    if flag == 1:
        outliers = find_outliers(B)
        B_stationary[outliers] = np.nan
        B_stationary_no_outliers = np.isnan(B_stationary).astype(int)
        no_outlier_sequences = zero_runs(B_stationary_no_outliers)

        ## remove sequences that are too short
        lst = no_outlier_sequences
        temp_list = [x for x in lst if x[1] - x[0] >= 30]
        no_outlier_sequences = temp_list
        n_outliers = np.zeros(shape=(20,))
        n_outliers[1] = 1
        count = 0
        count_iter = 0
        while len(outliers) > 0:
            #print count_iter
            if len(np.unique(n_outliers)) == 1 or count_iter>=max_iter:
                break
            if count < len(n_outliers):
                n_outliers[count] = len(outliers)
                count += 1
            else:
                count = 0
            for o in outliers:
                predicted = 0
                forward_seq = [x for x in no_outlier_sequences if x[1] <= o]
                if forward_seq:
                    temp_stationary = B_stationary.copy()
                    forward_seq = forward_seq[-1]  # pick the closest sequence to the outlier
                    forward_seq_data = B_stationary[forward_seq[0]:forward_seq[1]]
                    try:
                        AR_f, p, f = arburg(forward_seq_data, order=len(forward_seq_data) - 1, criteria='AICc')
                    except:
                        break
                    m = len(AR_f)
                    predicted_f = 0.0
                    for i in range(forward_seq[1], o + 1):
                        for j in range(m):
                            predicted_f += AR_f[j] * temp_stationary[i - 1 - j]
                        temp_stationary[i] = predicted_f
                    predicted_f = predicted_f + B_baseline[o]
                else:
                    predicted_f = 0

                backward_seq = [x for x in no_outlier_sequences if x[0] > o]
                if backward_seq:
                    temp_stationary = B_stationary.copy()
                    backward_seq = backward_seq[0]
                    backward_seq_data = B_stationary[backward_seq[0]:backward_seq[1]]
                    try:
                        AR_b, p, f = arburg(backward_seq_data, order=len(backward_seq_data) - 1, criteria='AICc')
                    except:
                        break
                    m = len(AR_b)
                    predicted_b = 0.0
                    for i in range(backward_seq[0] - 1, o - 1, -1):
                        for j in range(m):
                            predicted_b += AR_b[j] * temp_stationary[i + 1 + j]
                        temp_stationary[i] = predicted_b
                    predicted_b = predicted_b + B_baseline[o]
                else:
                    predicted_b = 0

                if predicted_b != 0 and predicted_f != 0:
                    predicted = (predicted_b + predicted_f) / 2
                elif predicted_b == 0:
                    predicted = predicted_f
                elif predicted_f == 0:
                    predicted = predicted_b
                else:
                    predicted_b = B[o]
                if np.abs(predicted) >= 200:  # and np.abs(predicted) > starting_B[o]:
                    B[o] = starting_B[o]
                elif not np.isnan(predicted):
                    B[o] = int(round(predicted))

            B_baseline = filtfilt(b, a, B)
            B_stationary = B - B_baseline
            outliers = find_outliers(B)
            B_stationary[outliers] = np.nan
            B_stationary_no_outliers = np.isnan(B_stationary).astype(int)
            no_outlier_sequences = zero_runs(B_stationary_no_outliers)

            ## remove sequences that are too short
            lst = no_outlier_sequences
            temp_list = [x for x in lst if x[1] - x[0] >= 30]
            no_outlier_sequences = temp_list
            count_iter += 1

    for i in range(len(B)):
        if B[i] < 140:
            B[i] = starting_B[i]
            flag = 1
        elif np.abs(B[i] - np.mean(starting_B)) > np.abs(starting_B[i] - np.mean(starting_B)):
            B[i] = starting_B[i]
    return B.astype(int)


def get_X_point(beat, B, C, fs=500, length=60, from_C=False, from_B=True):
    LVET_min = 100  # in samples
    LVET_max = (-1.3 * (60 * 2 * length / 1000) + 423) / 2  # in samples
    LVET_max2 = (-1.3 * (60 * 2 * length / 1000) + 483) / 2  # in samples
    d3zdt3_total = get_d3zdt3(beat, fs)
    d2_total = get_d2zdt2(beat, fs)
    if from_C:
        start = C + 25
        end = len(beat) - 15
    elif from_B:
        start = B + LVET_min
        end = B + LVET_max
        if end > len(beat) or end < start:
            end = end - 15
    start = int(start)
    end = int(end)
    start_2 = start
    end_2 = int(B + LVET_max2)
    B2X = beat[int(start):int(end)]
    if B2X.size and C < len(beat) - 5:
        if len(B2X) > 3:
            # look for all local min in B2X segment
            first_diff = np.ediff1d(B2X)
            local_min = argrelextrema(B2X, np.less)[0]
            d3zdt3 = d3zdt3_total[start:end]
            d3zdt3[-1] = 0
            for i in range(len(d3zdt3) - 1):
                if first_diff[i] < 0:
                    d3zdt3[i] = 0
            local_max = argrelextrema(d3zdt3, np.greater)[0]
            for i in range(len(local_max)):
                for j in range(15):
                    if local_max[i] + j < len(B2X):
                        if local_max[i] + j in local_min:
                            local_max[i] = local_max[i] + j
                            break
                    if local_max[i] - j in local_min:
                        local_max[i] = local_max[i] - j
                        break

            local_max = np.unique(local_max) + start
        else:
            d3zdt3 = d3zdt3_total[start_2:end_2]
            local_max = argrelextrema(d3zdt3, np.greater)[0] + start_2
        if not local_max.size:
            d3zdt3 = d3zdt3_total[start_2:end_2]
            local_max = argrelextrema(d3zdt3, np.greater)[0] + start_2
        if local_max.size > 1:
            local_max_height = beat[local_max]
            X = np.argmin(local_max_height)
            if local_max[X] > B + LVET_min:
                X = local_max[X]
            else:
                local_max = local_max[1:]
                local_max_height = beat[local_max]
                X = np.argmin(local_max_height)
                X = local_max[X]
        elif local_max.size == 1:
            X = local_max
        else:
            error = 1
            X = end - 1
    else:
        X = start + 1
        error = 1

    return X
    # in final version raise and track error


def get_correct_X(X, B, length):
    import matplotlib.pyplot as plt
    starting_X = X.copy()
    min_X = B + 100
    # max_X = B + (-1.3*(60*2*length/1000) + 393)/2
    max_X = B + (-1.3 * (60 * 2 * length / 1000) + 473) / 2
    b, a = butter(4, 0.1)
    X_baseline = filtfilt(b, a, X)
    X_stationary = X - X_baseline

    ## based on outliers break the series into chunks with no outliers
    outliers = find_outliers(X)
    X_stationary[outliers] = np.nan
    X_stationary_no_outliers = np.isnan(X_stationary).astype(int)
    no_outlier_sequences = zero_runs(X_stationary_no_outliers)

    ## remove sequences that are too short
    lst = no_outlier_sequences
    temp_list = [x for x in lst if x[1] - x[0] >= 30]
    no_outlier_sequences = temp_list
    n_outliers = np.zeros(shape=(20,))
    n_outliers[1] = 1
    count = 0
    count_iter = 0
    max_iter = 20
    while len(outliers) > 0:
        #print count_iter
        if len(np.unique(n_outliers)) == 1 or count_iter >= max_iter:
            break
        if count < len(n_outliers):
            n_outliers[count] = len(outliers)
            count += 1
        else:
            count = 0
        for o in outliers:
            predicted = 0
            forward_seq = [x for x in no_outlier_sequences if x[1] <= o]
            if forward_seq:
                temp_stationary = X_stationary.copy()
                forward_seq = forward_seq[-1]  # pick the closest sequence to the outlier
                forward_seq_data = X_stationary[forward_seq[0]:forward_seq[1]]
                try:
                    AR_f, p, f = arburg(forward_seq_data, order=len(forward_seq_data) - 1, criteria='AICc')
                except:
                    break
                m = len(AR_f)
                predicted_f = 0.0
                for i in range(forward_seq[1], o + 1):
                    for j in range(m):
                        predicted_f += AR_f[j] * temp_stationary[i - 1 - j]
                    temp_stationary[i] = predicted_f
                predicted_f = predicted_f + X_baseline[o]
            else:
                predicted_f = 0

            backward_seq = [x for x in no_outlier_sequences if x[0] > o]
            if backward_seq:
                temp_stationary = X_stationary.copy()
                backward_seq = backward_seq[0]
                backward_seq_data = X_stationary[backward_seq[0]:backward_seq[1]]
                try:
                    AR_b, p, f = arburg(backward_seq_data, order=len(backward_seq_data) - 1, criteria='AICc')
                except:
                    break
                m = len(AR_b)
                predicted_b = 0.0
                for i in range(backward_seq[0] - 1, o - 1, -1):
                    for j in range(m):
                        predicted_b += AR_b[j] * temp_stationary[i + 1 + j]
                    temp_stationary[i] = predicted_b
                predicted_b = predicted_b + X_baseline[o]
            else:
                predicted_b = 0

            if predicted_b != 0 and predicted_f != 0:
                predicted = (predicted_b + predicted_f) / 2
            elif predicted_b == 0:
                predicted = predicted_f
            elif predicted_f == 0:
                predicted = predicted_b
            # else:
                #print ('clean up in aisle outlier!!')
            if np.abs(predicted) >= max_X[o]:  # and np.abs(predicted) > starting_B[o]:
                X[o] = starting_X[o]
            else:
                X[o] = int(round(predicted))

        X_baseline = filtfilt(b, a, X)
        X_stationary = X - X_baseline
        outliers = find_outliers(X)
        X_stationary[outliers] = np.nan
        X_stationary_no_outliers = np.isnan(X_stationary).astype(int)
        no_outlier_sequences = zero_runs(X_stationary_no_outliers)

        ## remove sequences that are too short
        lst = no_outlier_sequences
        temp_list = [x for x in lst if x[1] - x[0] >= 30]
        no_outlier_sequences = temp_list
        count_iter +=1
    flag = 0
    for i in range(len(X)):
        if X[i] < min_X[i]:
            X[i] = starting_X[i]
            flag = 1
        elif np.abs(X[i] - np.mean(starting_X)) > np.abs(starting_X[i] - np.mean(starting_X)):
            X[i] = starting_X[i]
            flag = 1

    if flag == 1:
        outliers = find_outliers(X)
        X_stationary[outliers] = np.nan
        B_stationary_no_outliers = np.isnan(X_stationary).astype(int)
        no_outlier_sequences = zero_runs(X_stationary_no_outliers)

        ## remove sequences that are too short
        lst = no_outlier_sequences
        temp_list = [x for x in lst if x[1] - x[0] >= 30]
        no_outlier_sequences = temp_list
        n_outliers = np.zeros(shape=(20,))
        n_outliers[1] = 1
        count = 0
        count_iter = 0
        while len(outliers) > 0:
            #print count_iter
            if len(np.unique(n_outliers)) == 1 or count_iter >= max_iter:
                break
            if count < len(n_outliers):
                n_outliers[count] = len(outliers)
                count += 1
            else:
                count = 0
            for o in outliers:
                predicted = 0
                forward_seq = [x for x in no_outlier_sequences if x[1] <= o]
                if forward_seq:
                    temp_stationary = X_stationary.copy()
                    forward_seq = forward_seq[-1]  # pick the closest sequence to the outlier
                    forward_seq_data = X_stationary[forward_seq[0]:forward_seq[1]]
                    try:
                        AR_f, p, f = arburg(forward_seq_data, order=len(forward_seq_data) - 1, criteria='AICc')
                    except:
                        break
                    m = len(AR_f)
                    predicted_f = 0.0
                    for i in range(forward_seq[1], o + 1):
                        for j in range(m):
                            predicted_f += AR_f[j] * temp_stationary[i - 1 - j]
                        temp_stationary[i] = predicted_f
                    predicted_f = predicted_f + X_baseline[o]
                else:
                    predicted_f = 0

                backward_seq = [x for x in no_outlier_sequences if x[0] > o]
                if backward_seq:
                    temp_stationary = X_stationary.copy()
                    backward_seq = backward_seq[0]
                    backward_seq_data = X_stationary[backward_seq[0]:backward_seq[1]]
                    try:
                        AR_b, p, f = arburg(backward_seq_data, order=len(backward_seq_data) - 1, criteria='AICc')
                    except:
                        break
                    m = len(AR_b)
                    predicted_b = 0.0
                    for i in range(backward_seq[0] - 1, o - 1, -1):
                        for j in range(m):
                            predicted_b += AR_b[j] * temp_stationary[i + 1 + j]
                        temp_stationary[i] = predicted_b
                    predicted_b = predicted_b + X_baseline[o]
                else:
                    predicted_b = 0

                if predicted_b != 0 and predicted_f != 0:
                    predicted = (predicted_b + predicted_f) / 2
                elif predicted_b == 0:
                    predicted = predicted_f
                elif predicted_f == 0:
                    predicted = predicted_b
                else:
                    predicted_b = X[o]
                if np.abs(predicted) >= max_X[o]:  # and np.abs(predicted) > starting_B[o]:
                    X[o] = starting_X[o]
                else:
                    X[o] = int(round(predicted))

            B_baseline = filtfilt(b, a, X)
            B_stationary = X - B_baseline
            outliers = find_outliers(X)
            B_stationary[outliers] = np.nan
            B_stationary_no_outliers = np.isnan(B_stationary).astype(int)
            no_outlier_sequences = zero_runs(B_stationary_no_outliers)

            ## remove sequences that are too short
            lst = no_outlier_sequences
            temp_list = [x for x in lst if x[1] - x[0] >= 30]
            no_outlier_sequences = temp_list
            count_iter+=1

    for i in range(len(X)):
        if X[i] < min_X[i]:
            X[i] = starting_X[i]
            flag = 1
        elif np.abs(X[i] - np.mean(starting_X)) > np.abs(starting_X[i] - np.mean(starting_X)):
            X[i] = starting_X[i]

    return X.astype(int)

