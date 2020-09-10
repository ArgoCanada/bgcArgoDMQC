#!/usr/bin/python

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from bgcArgo import sprof, profiles

sns.set(context='paper', style='whitegrid', palette='colorblind')
sage_path = Path('/Users/gordonc/Documents/data/Argo/sage/')
sage_files = list(sage_path.glob('*gui*.h5'))

for fn in sage_files:
    wmo = int(str(fn).split('\\')[-1].split('_')[-1].split('.')[0])

    datapath = Path('/Users/gordonc/Documents/data/')
    profiles.set_dirs(argo_path=datapath / 'Argo', woa_path=datapath / 'WOA18', ncep_path=datapath / 'NCEP')
    syn = sprof(wmo)
    gains = syn.calc_gains(ref='WOA')

    sf = pd.DataFrame(dict(pyCYCLE=syn.CYCLE, pyGAINS=gains))
    mf = pd.read_hdf(fn)
    
    df = pd.merge(sf, mf, left_on='pyCYCLE', right_on='CYCLE')
    print(df.head())
