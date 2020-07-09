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
    ### do save steps ###

    # load data
    data = loadmat(gfn)

    # extract relevant data
    track = data['track']
    cycle = np.squeeze(data['cycles'])
    gains = np.squeeze(data['GAINS'])
    surf_sat = data['SURF_SAT'][:,1]
    woa_sat  = data['refdata'][:,3]

    # specific variables
    time = track[:,0]
    lon = track[:,2]
    lat = track[:,3]

    # dict for dateframe
    df_dict = dict(SDN=time, CYCLE=cycle, LAT=lat, LON=lon, SURF_SAT=surf_sat, REF_SAT=woa_sat, GAINS=gains)
    # dataframe 
    df = pd.DataFrame(df_dict)
