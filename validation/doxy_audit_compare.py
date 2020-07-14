#!/usr/bin/python

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import bgcArgo as bgc

# load in DOXY audit file - most recent one was on July 9 2020
audit_file = Path('../data/DOXY_audit_070720.TXT')
df = pd.read_csv(audit_file, sep='\t', header=25)

# download the synthetic, meta, and BRtraj files for each float in the audit
dacpath = '/ifremer/argo/dac'
fltpath = ['{}/{}/{}'.format(dacpath, bgc.get_dac(w), w) for w in df.WMO.unique()]
local_path = '/Users/gordonc/Documents/data/Argo'
# bgc.io.get_argo(fltpath, local_path=local_path, mode='summary')
bgc.set_dirs(argo_path='/Users/gordonc/Documents/data/Argo', woa_path='/Users/gordonc/Documents/data/WOA18')

for wmo in df.WMO.unique():
    sub = df[df.WMO == wmo]
    syn = bgc.sprof(wmo)
    syn.clean()
    syn.calc_gains(ref='WOA')
    for i in range(sub.shape[0]):
        cycle = sub.cycle.iloc[i]
        ix = syn.CYCLE == cycle
        if any(ix):
            print('gain from python bgcArgo: {:2f}, gain from matlab SAGE-O2: {:2f}'.format(syn.gains[ix], sub.G_raw.iloc[i]))
        else:
            print('Cycle number not present in synthetic file')
