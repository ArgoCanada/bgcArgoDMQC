#!/usr/bin/python

from pathlib import Path

import numpy as np
import pandas as pd

from bgcArgo import sprof

# setup directories
datadir   = Path('/Users/gordonc/Documents/data')
woa_path  = datadir / 'WOA18'
ncep_path = datadir / 'NCEP'
argo_path = datadir / 'Argo/meds'
sprof.set_dirs(woa_path=woa_path, ncep_path=ncep_path, argo_path=argo_path)

# create synthetic profile object
flt = sprof(4902481)

# run through all required data
flt.get_track()
flt.get_woa()
flt.get_ncep()

# create dataframe for all relevant data
df = pd.DataFrame()
# track
df['time'] = flt.track[:,0]
df['lat']  = flt.track[:,1]
df['lon']  = flt.track[:,2]
# reference data
df['ncep'] = flt.NCEP
df['woa']  = flt.WOA[0,:]
# WOA reference weights
df['woa_time_wt'] = flt.__WOAweights__[0]
df['woa_lat_wt']  = flt.__WOAweights__[1]
df['woa_lon_wt']  = flt.__WOAweights__[2]

# NCEP weights
nf = pd.DataFrame()
nf['ncep_time_wt'] = flt.__NCEPweights__[0]
nf['ncep_lat_wt']  = flt.__NCEPweights__[1]
nf['ncep_lon_wt']  = flt.__NCEPweights__[2]

# save files to hdf
store = pd.HDFStore('../data/4902481_trackdata.h5')
store.put('df', df, data_columns=df.columns)
store.close()

# save files to hdf
store = pd.HDFStore('../data/4902481_ncep_weights.h5')
store.put('df', nf, data_columns=df.columns)
store.close()
