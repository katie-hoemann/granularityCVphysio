from __future__ import division
import numpy as np

def pep_max_min_checks(pep,max=200,min=30):
    checks = np.zeros(shape=len(pep))
    for i in range(len(pep)):
        if pep[i] >=max or pep[i] <=min:
            checks[i] = 1
    return checks

def lvet_max_min_checks(lvet,max=500,min=100):
    checks = np.zeros(shape=len(lvet))
    for i in range(len(lvet)):
        if lvet[i] >=max or lvet[i] <=min:
            checks[i] = 1
    return checks
