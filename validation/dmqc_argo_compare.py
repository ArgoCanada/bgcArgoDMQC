#!/usr/bin/python

import sys
from pathlib import Path

from netCDF4 import Dataset

import numpy as np
import pandas as pd
import bgcArgo as bgc

argopath = '/Users/gordonc/Documents/data/Argo/'
bgc.set_dirs(argo_path=argopath, woa_path='/Users/gordonc/Documents/data/WOA18/')
wmos = [3900345]

for wmo in wmos:
    syn = bgc.sprof(wmo)
    gains = syn.calc_gains(ref='WOA')

    files = bgc.get_files(local_path=argopath, wmo_numbers=[wmo], mission='B', mode='D')

    file_time  = np.array(len(files)*[np.nan])
    file_gains = np.array(len(files)*[np.nan])
    file_msgs  = np.array(len(files)*[256*' '])
    syn_gains  = np.array(len(files)*[np.nan])

    for i,fn in enumerate(files):
        nc = Dataset(Path(fn))
        gain, msg = bgc.util.read_gain_value(nc)
        c = nc.variables['CYCLE_NUMBER'][:][0]
        if c > syn.CYCLE[-1]:
            syn_gain = np.nan
        else:
            syn_gain = gains[syn.CYCLE == c][0]

        file_time[i]  = nc.variables['JULD'][:][0]
        file_gains[i] = gain
        file_msgs[i]  = msg
        syn_gains[i]  = syn_gain
    
    diffs = syn_gains - file_gains
    