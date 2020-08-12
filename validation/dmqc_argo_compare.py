#!/usr/bin/python

import sys
from pathlib import Path

from netCDF4 import Dataset

import numpy as np
import pandas as pd
import bgcArgo as bgc

def get_all_wmos():

    argo_path = Path('/Users/gordonc/Documents/data/Argo/')
    wmos = []
    dacs = list(argo_path.glob('*'))
    dacs.pop(-1)
    for dac in dacs:
        wmos = wmos + list(dac.glob('*'))

    wmos = [int(p.as_posix().split('/')[-1]) for p in wmos]
    
    return wmos

argopath = '/Users/gordonc/Documents/data/Argo/'
bgc.set_dirs(argo_path=argopath, woa_path='/Users/gordonc/Documents/data/WOA18/')
wmos = get_all_wmos()

for wmo in wmos:
    files = bgc.get_files(local_path=argopath, wmo_numbers=[wmo], mission='B', mode='D')

    if len(files) > 0:
        syn = bgc.sprof(wmo)
        gains = syn.calc_gains(ref='WOA')

        file_time  = np.array(len(files)*[np.nan])
        # file_gains = np.array(len(files)*[np.nan])
        # file_msgs  = np.array(len(files)*[256*' '])
        file_gains = []
        file_msgs  = []
        syn_gains  = np.array(len(files)*[np.nan])

        for i,fn in enumerate(files):
            nc = Dataset(Path(fn))
            gain, msg = bgc.util.read_gain_value(nc)
            print(gain)
            print(msg)
            c = nc.variables['CYCLE_NUMBER'][:][0]
            if c > syn.CYCLE[-1]:
                syn_gain = np.nan
            else:
                syn_gain = gains[syn.CYCLE == c][0]

            file_time[i]  = nc.variables['JULD'][:][0]
            file_gains.append(gain)
            file_msgs.append(msg)
            syn_gains[i]  = syn_gain
        
        # diffs = syn_gains - file_gains
