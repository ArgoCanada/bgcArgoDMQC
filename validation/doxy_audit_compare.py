#!/usr/bin/python

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import bgcArgoDMQC as bgc

# load in DOXY audit file - most recent one was on July 9 2020
audit_file = Path('/Users/gordonc/Documents/argo/doxy-audit/DOXY_audit_112020.TXT')
df = pd.read_csv(audit_file, sep='\t', header=25)

# download the synthetic, meta, and BRtraj files for each float in the audit
local_path = '/Users/gordonc/Documents/data/Argo'
# for wmo in df.WMO.unique():
    # bgc.io.get_argo(wmo, local_path=local_path, ftype='summary')

with open(Path('../data/doxy_audit_vs_bgcArgo_py_comparison_20201230.csv'),'w') as fid:
    fid.write('WMO,CYCLE,DAC,DATE,pyGAIN,sageGAIN')
    for wmo in df.WMO.unique():
        sub = df[df.WMO == wmo]
        syn = bgc.sprof(wmo)
        syn.calc_gains(ref='WOA')
        for i in range(sub.shape[0]):
            cycle = sub.cycle.iloc[i]
            ix = syn.CYCLE == cycle
            print('(pyGain, matlabGain) = ({:.2f}, {:.2f})'.format(syn.gains[ix][0], sub.G_raw.iloc[i]))
            fid.write('\n{:d},{:d},{:s},{:s},{:.2f},{:.2f}'.format(wmo,cycle,sub.DAC.iloc[i],sub['profile date'].iloc[i],syn.gains[ix][0],sub.G_raw.iloc[i]))
