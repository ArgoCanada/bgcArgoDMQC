#!/usr/bin/python

from pathlib import Path

from scipy.io import loadmat

import numpy as np
import pandas as pd

# get sage-produced matlab float files
matlab_file_path = '/Users/gordonc/Documents/data/Argo/sage/'
float_files = list(Path(matlab_file_path).glob('float*.mat'))
gui_files = list(Path(matlab_file_path).glob('gui*.mat'))

# loop through files and save relevant info to dataframes
for ffn, gfn in zip(float_files, gui_files):
    # load data
    data = loadmat(ffn)
    fltID = str(ffn).split('\\')[-1].split('.')[0].split('_')[-1]

    # extract data from matlab dict
    flttype = data['floatTYPE'][0]
    cycles  = np.squeeze(data['cycles'])
    track   = data['track']
    pts     = data['PTSdata']
    oxy     = data['O2data'][0,0]

    # specific variables from sprof
    time = oxy[:,0]
    cycle = oxy[:,1]
    PSAL = oxy[:,2]
    PRES = oxy[:,3]
    TEMP = oxy[:,4]
    DOXY = oxy[:,6]
    O2SAT = oxy[:,10]

    # dict for dataframe
    df_dict = dict(SDN=time, CYCLE=cycle, PRES=PRES, PSAL=PSAL, TEMP=TEMP, DOXY=DOXY, O2SAT=O2SAT)
    # dataframe
    df = pd.DataFrame(df_dict)

    # save dataframe
    store = pd.HDFStore(Path('/Users/gordonc/Documents/data/Argo/sage/df_floatdata_{}.h5'.format(fltID)))
    store.put('df', df, data_columns=df.columns)
    store.close()

    # load data
    data = loadmat(gfn)
    fltID = str(gfn).split('\\')[-1].split('.')[0].split('_')[-1]

    # extract relevant data
    track = data['track']
    cycle = np.squeeze(data['cycles'])
    gains = np.squeeze(data['GAINS'])
    surf_sat = data['SURF_SAT'][:,1]
    woa_sat  = data['refdata'][:,3]

    if gains.shape[0] != surf_sat.shape[0]:
        print('gain shape ({}) does not match data shape ({}) for float {}, calculating independently using data'.format(gains.shape[0], surf_sat.shape[0], fltID))
        gains = woa_sat/surf_sat

    # specific variables
    time = track[:,0]
    lon = track[:,2]
    lat = track[:,3]

    # dict for dateframe
    df_dict = dict(SDN=time, CYCLE=cycle, LAT=lat, LON=lon, SURF_SAT=surf_sat, REF_SAT=woa_sat, GAINS=gains)
    # dataframe 
    df = pd.DataFrame(df_dict)

    # save dataframe
    store = pd.HDFStore(Path('/Users/gordonc/Documents/data/Argo/sage/df_guidata_{}.h5'.format(fltID)))
    store.put('df', df, data_columns=df.columns)
    store.close()
