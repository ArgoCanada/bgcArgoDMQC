#!/usr/bin/python

from pathlib import Path

import numpy as np
import pandas as pd

from bgcArgo import sprof
sprof.set_dirs(
    argo_path='/Users/gordonc/Documents/data/Argo', 
    woa_path='/Users/gordonc/Documents/data/WOA18',
    ncep_path='/Users/gordonc/Documents/data/NCEP'
)

sagepath = Path('/Users/gordonc/Documents/data/Argo/sage/')
gui_files = sagepath.glob('*guidata*.h5')
flt_files = sagepath.glob('*floatdata*.h5')

df = pd.DataFrame(columns=['SDN', 'CYCLE', 'LAT', 'LON', 'SURF_SAT', 'REF_SAT', 'GAINS', 'pyGAIN', 'WMO'])

for gfn, ffn in zip(gui_files, flt_files):
    gdf = pd.read_hdf(gfn)
    fdf = pd.read_hdf(ffn)

    wmo = int(str(ffn).split('\\')[-1].split('_')[-1].split('.')[0])
    syn = sprof(wmo)
    woa_gains = syn.calc_gains(ref='WOA')
    try:
        print(wmo)
        ncep_gains = syn.calc_gains()
        ncep_flag = True
    except:
        print('shoot {}'.format(wmo))
        ncep_flag = False

    pyWOAGAIN = np.nan*np.ones((gdf.shape[0]))
    pyAIRGAIN = np.nan*np.ones((gdf.shape[0]))

    for i,c in enumerate(gdf.CYCLE):
        if c in syn.CYCLE:
            ix = syn.CYCLE == c
            g = woa_gains[ix]
            if ncep_flag:
                a = ncep_gains[ix]
            if g.shape[0] > 1:
                g = np.nanmean(g)
                if ncep_flag:
                    a = np.nanmean(a)
            else:
                g = g[0]
                if ncep_flag:
                    a = a[0]
            pyWOAGAIN[i] = g
            if ncep_flag:
                pyAIRGAIN[i] = a

    gdf['pyWOAGAIN'] = pyWOAGAIN
    gdf['pyAIRGAIN'] = pyAIRGAIN
    gdf['WMO'] = np.array(gdf.shape[0]*[wmo])

    df = pd.concat([df, gdf])

store = pd.HDFStore(Path('../data/sage_gui_comparison.h5'))
store.put('df', df, data_columns=df.columns)
store.close()