import numpy as np
import os
def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]


def zero_runs(a):
    iszero = np.concatenate(([0], np.equal(a, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges

def extract_numbers_from_strings(input):
    import re
    return re.findall(r'\d+', str(input))

def process_times(sitting_times,prompt_time=0,start=300,end=30):

    prompt_time = 0
    start = prompt_time + start
    end = prompt_time - end
    total_sitting = 0
    if not sitting_times.size:
        total_sitting = 0
    else:
        if len(sitting_times.shape) > 1:
            for times in sitting_times:
                if times[0]<=start and times[1]>=end:
                    total_sitting = total_sitting + times[0] - times[1]
                elif times[0]>start and start>=times[1]>=end:
                    total_sitting = total_sitting + start - times[1]
                elif end<=times[0]<=start and times[1]<end:
                    total_sitting = total_sitting + times[0] - end
                elif times[0]>start and times[1]<end:
                    total_sitting = total_sitting + start - end
                    break
                elif times[0]>start and times[1]>start:
                    total_sitting = total_sitting
                elif times[0]<end and times[1]<end:
                    total_sitting = total_sitting
            total_sitting = 100 * total_sitting / (start - end)
        else:
            times = sitting_times
            if times[0] <= start and times[1] >= end:
                total_sitting = total_sitting + times[0] - times[1]
            elif times[0] > start and start >= times[1] >= end:
                total_sitting = total_sitting + start - times[1]
            elif end <= times[0] <= start and times[1] < end:
                total_sitting = total_sitting + times[0] - end
            elif times[0] > start and times[1] < end:
                total_sitting = total_sitting + start - end
            elif times[0] > start and times[1] > start:
                total_sitting = total_sitting
            elif times[0] < end and times[1] < end:
                total_sitting = total_sitting
            total_sitting = 100 * total_sitting / (start - end)
            # if total_sitting!=0 and total_sitting!=100:
            #     total_sitting = 100 - total_sitting
        # if total_sitting == 0:
        #     total_sitting = 100

    return np.round(total_sitting,2)

def get_immediate_subdirectories(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]
