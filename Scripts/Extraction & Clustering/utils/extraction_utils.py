from __future__ import division
import numpy as np
import os
import re

def convert_to_seconds(time,twelve_hour=False):
    """
    this function converts time in 12/24 hours format to seconds from midnight
    :param time: input time in twelve 12/24 hour format
    :param twelve_hour: True if 12 hour format is input
    :return: ms/1000: Time in seconds from midnight
    """
    if twelve_hour:
        if time[-2] == 'A':
            hh = int(time[:2])
            if hh == 12:
                hh = 0
            mm = int(time[3:5])
            ss = int(time[6:8])
            ms = int(time[9:12])
            ms = ms + hh*60*60*1000 + mm*60*1000 + ss*1000
        else:
            hh = int(time[:2])
            if hh!=12:
                hh = hh + 12
            mm = int(time[3:5])
            ss = int(time[6:8])
            ms = int(time[9:12])
            ms = ms + hh*60*60*1000 + mm*60*1000 + ss*1000
    else:
        hh = int(time[:2])
        mm = int(time[3:5])
        ss = int(time[6:8])
        us = int(time[9:])
        ms = us/1000
        ms = ms + hh*60*60*1000 + mm*60*1000 + ss*1000
    return ms/1000


def load_event(event_file_name,day,subject_id,drop_directory,threshold=30,count=0):
    """
    This function takes the physio information for each sitting segment and puts it in a separate file
    :param event_file_name: Mindware Event File
    :param day: current day being processed
    :param subject_id: current participant being processed
    :param drop_directory: drop extracted raw data in the same column ordering as Mindware Physio files
    :param threshold: ignore events SHORTER than this threshold (in seconds)
    :return: sitting_lengths: length of sitting segments in current event file, also write all sitting segments
    longer than the threshold to file
    """
    f = open(event_file_name, 'r')
    l = f.readlines()
    f.close()
    acquisition_start_time = l[1].split('\t')[3][:-1]
    acquisition_start_date = l[1].split('\t')[2][3:5]
    acquisition_start_time = convert_to_seconds(acquisition_start_time,twelve_hour=True)

    sitting_start_times = []
    sitting_end_times = []
    for (itr,m) in enumerate(l):
        if len(m.split(' => ')) > 1:
            if m.split(' => ')[1][:7] == 'Sitting':
                sitting_temp = (m.split('\t')[3].replace('\n',''),m.split('\t')[2].replace('\n',''))
                # for i in range(itr+1,len(l)):
                if itr != len(l)-1:
                    if l[itr+1].split('\t')[1][0] == 'P':
                        if l[itr+1].split(' => ')[0][-7:] == 'Sitting':
                            sitting_start_times.append(sitting_temp)
                            sitting_end_times.append((l[itr+1].split('\t')[3].replace('\n',''),l[itr+1].split('\t')[2].replace('\n','')))
                            # break
    sitting_length = []
    for (start,end) in zip(sitting_start_times,sitting_end_times):
        start_time = start[0]
        start_date = start[1]
        end_time = end[0]
        end_date = end[1]
        if start_date == end_date:
            current_length = convert_to_seconds(end_time,twelve_hour=True) \
                             - convert_to_seconds(start_time,twelve_hour=True)
        else:
            current_length = convert_to_seconds(end_time, twelve_hour=True) \
                             - convert_to_seconds(start_time, twelve_hour=True) + 86400
        if current_length >= threshold:
            count += 1
            start_time = convert_to_seconds(start_time,twelve_hour=True)
            end_time = convert_to_seconds(end_time,twelve_hour=True)
            sitting_length.append(current_length)
            if start_date[3:5] == acquisition_start_date:
                start_trigger = start_time -  acquisition_start_time
            else:
                start_trigger = start_time - acquisition_start_time + (int(start_date[3:5]) - int(acquisition_start_date)) * 86400
            start_trigger = int(start_trigger * 1000)
            if start_trigger % 2!=0:
                start_trigger = start_trigger + 1
                start_trigger = start_trigger/1000
            else:
                start_trigger = start_trigger/1000

            if end_date[3:5] == acquisition_start_date:
                end_trigger = end_time -  acquisition_start_time
            else:
                end_trigger = end_time - acquisition_start_time + (int(end_date[3:5]) - int(acquisition_start_date)) * 86400
            end_trigger = int(end_trigger * 1000)
            if end_trigger % 2!=0:
                end_trigger = end_trigger + 1
                end_trigger = end_trigger/1000
            else:
                end_trigger = end_trigger/1000
            physio_file = event_file_name.replace("_event.txt",".txt")
            export_event(physio_file=physio_file,starting_time=start_trigger,
                         ending_time=end_trigger,day=day,subject_id=subject_id,event_id=count,path=drop_directory)
    return sitting_length

def seek_to_position(fp, starting_time):
    """
    This function finds the starting point of an event in the big physio text files
    :param fp: file pointer
    :param starting_time: time the event starts (in seconds, as recorded in the mindware physio file first column)
    :return:
    """
    x1 = 1000
    x2 = 2000

    fp.seek(0, 2)
    last_pointer = fp.tell()

    fp.seek(x1)
    fp.readline()
    y1 = float(fp.readline().split('\t')[0])

    count = 0
    while abs(starting_time - y1 - 0.0005) > 0.001 and count < 30:
        count += 1
        fp.seek(x2)
        fp.readline()
        y2 = float(fp.readline().split('\t')[0])


        x1 = x1 + int((starting_time - y1) * (x2 - x1) / (y2 - y1))
        if x1 > last_pointer:
            fp.seek(x1)
            break

        fp.seek(x1)
        fp.readline()
        test = fp.readline()
        if test == '':
            fp.seek(last_pointer + 100)  # there's an error
            return

        fp.seek(x1)
        fp.readline()
        kkk = fp.readline()

        y1 = float(kkk.split('\t')[0])
        x2 = x1 + 1000


def export_event(physio_file,
                 starting_time,
                 ending_time,
                 day = 0,
                 subject_id = None,
                 event_id = 1,
                 path=None,):
    """

    :param physio_file: mindware physio file to extract event from
    :param starting_time: event starting time
    :param ending_time: event ending time
    :param day: current day
    :param subject_id: current subject
    :param event_id: event counter for the day
    :param path: parent directory for dumping the extracted data
    :return:
    """
    events_directory = path + subject_id + '/' + day + '/' + 'events/'
    if not os.path.exists(events_directory):
        os.makedirs(events_directory)
    event_f = open(events_directory + str(subject_id) + '_day_' +
                   str(np.array(re.findall('\d+',day[3:6])).astype(int)[0]) +
                   '_event_' + str(event_id) + '_tick_' + str(starting_time).replace(".", "_") + '.txt',
                   'w')
    print (event_f)
    fp = open(physio_file, 'r')
    fp.seek(0, 2)
    EOF_pointer = fp.tell()
    seek_to_position(fp, ending_time)
    end_pointer = fp.tell()
    seek_to_position(fp, starting_time)
    nxLine = 'start'
    new_time = starting_time
    temp_ECG = []
    temp_IMP = []
    temp_EDA = []
    while new_time < ending_time and nxLine != '':
        nxLine = fp.readline()
        if nxLine == '':
            break
        else:
            ln_split = nxLine.strip().split('\t')
            new_time = float(ln_split[0])
            temp_ECG.append(float(ln_split[1]))
            temp_IMP.append(float(ln_split[3]))
            temp_EDA.append(float(ln_split[4]))
            event_f.write(str(ln_split[0]) + '\t' + str(ln_split[1]) + '\t' + str(ln_split[2]) + '\t'
                          + str(ln_split[3]) + '\t' + str(ln_split[4]) + '\t' + str(ln_split[5]) + '\t'
                          + str(ln_split[6]) + '\t' + str(ln_split[7]) + '\n')

    fp.close()

