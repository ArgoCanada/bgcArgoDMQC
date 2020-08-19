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

def decode_file_gain_strings(fg_arr, msg_arr):

    GAIN = np.array(fg_arr.shape[0]*[np.nan])
    MSGS = np.array(fg_arr.shape[0]*[256*' '])

    for n,glist in enumerate(fg_arr):
        if type(glist) is float:
            g = glist
        else:
            if len(glist) == 1:
                if type(glist[0]) is str:
                    g = float(glist[0].split('=')[-1].strip())
                else:
                    g = np.nan
                msg = msg_arr[n][0]
            else:
                g = len(glist) * [np.nan]
                for i,gstr in enumerate(glist):
                    g[i] = float(gstr.split('=')[-1].strip())
                g = np.array(g)
                if (g == g[0]).all():
                    g = g[0]
                    msgs = msg_arr[n][0]
                else:
                    print('Different calibration values, taking last one')
                    g = g[-1]
                    msg = msg_arr[n][-1]
        
        GAIN[n] = g
        MSGS[n] = msg
    
    return GAIN, MSGS

argopath = '/Users/gordonc/Documents/data/Argo/'
bgc.set_dirs(argo_path=argopath, woa_path='/Users/gordonc/Documents/data/WOA18/')
wmos = get_all_wmos()

index = bgc.get_index()
wmos = pd.DataFrame(dict(wmos=wmos))
wmos = wmos[wmos.wmos.isin(index.wmo)].wmos.unique()

file_gains = []
file_msgs  = []
file_time  = np.array([])
syn_gains  = np.array([])
mean_gain  = np.array([])
arr_sdn    = np.array([])
arr_dac    = np.array([])
arr_wmo    = np.array([], dtype=int)
arr_cyc    = np.array([])

for wmo in wmos:
    print(wmo)
    files = bgc.get_files(local_path=argopath, wmo_numbers=[wmo], mission='B', mode='D')
    if len(files) > 0:
        syn = bgc.sprof(wmo)
        syn.clean()
        gains = syn.calc_gains(ref='WOA')

        sub_file_time  = np.array(len(files)*[np.nan])
        # sub_file_gains = np.array(len(files)*[np.nan])
        # sub_file_msgs  = np.array(len(files)*[256*' '])
        sub_syn_gains  = np.array(len(files)*[np.nan])
        sub_cycles     = np.array(len(files)*[np.nan])
        sub_wmo        = np.array(len(files)*[int(wmo)], dtype=int)
        sub_dac        = np.array(len(files)*[bgc.get_dac(wmo)])
        sub_sdn        = np.array(len(files)*[np.nan])

        for i,fn in enumerate(files):
            nc = Dataset(Path(fn))
            gain, msg = bgc.util.read_gain_value(nc)
            c = nc.variables['CYCLE_NUMBER'][:][0]
            sub_cycles[i] = c
            if c not in syn.CYCLE:
                syn_gain = np.nan
                syn_sdn  = np.nan
            else:
                syn_gain = gains[syn.CYCLE == c][0]
                syn_sdn  = syn.SDN[syn.CYCLE == c][0]

            sub_file_time[i]  = nc.variables['JULD'][:][0]
            file_gains.append(gain)
            file_msgs.append(msg)
            sub_syn_gains[i]  = syn_gain
            sub_sdn[i] = syn_sdn

        file_time = np.append(file_time, sub_file_time)
        syn_gains = np.append(syn_gains, sub_syn_gains)
        mean_gain = np.append(mean_gain, np.array(sub_syn_gains.shape[0]*[np.nanmean(sub_syn_gains)]))
        arr_wmo   = np.append(arr_wmo, sub_wmo)
        arr_dac   = np.append(arr_dac, sub_dac)
        arr_sdn   = np.append(arr_sdn, sub_sdn)
        arr_cyc   = np.append(arr_cyc, sub_cycles)
        
file_gains = np.array(file_gains, dtype=object)
file_msgs  = np.array(file_msgs, dtype=object)

processed_gains, processed_msgs = decode_file_gain_strings(file_gains, file_msgs)

df = pd.DataFrame(dict(WMO=arr_wmo, CYCLE=arr_cyc, DAC=arr_dac, DATE=arr_sdn, pyGAIN=syn_gains, pyMEANGAIN=mean_gain, argoGAIN=processed_gains, argoMSG=processed_msgs))

store = pd.HDFStore(Path('../data/argo_dmqc_comparison_20200819.h5'))
store.put('df', df, data_columns=df.columns)
store.close()
