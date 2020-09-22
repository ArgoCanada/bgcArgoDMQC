#!/usr/bin/python

import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from bgcArgo import sprof

plot_profiles = True

# summary comparison between bgcArgo and SAGE/DOXY audit
fn = Path('../data/doxy_audit_vs_bgcArgo_py_comparison_20200920.csv')
df = pd.read_csv(fn)
df['diffGAIN'] = np.abs(df.pyGAIN - df.sageGAIN)

audit_file = Path('../data/DOXY_audit_070720.TXT')
xf = pd.read_csv(audit_file, sep='\t', header=25)

nan  = df[df.pyGAIN.isna()]
big  = df[df.diffGAIN >= 0.2]
good = df[df.diffGAIN == 0.0]

Nnan = nan.shape[0]
if Nnan == 0:
    Nnan = 1

# by dac
sys.stdout.write('%nan (N)\t%big (N)\t%audit (N)\t dac\n')
for dac in xf.DAC.unique():
    sub_nan = nan[nan.DAC == dac]
    sub_big = big[big.DAC == dac]
    sub_audit = xf[xf.DAC == dac]
    sys.stdout.write('{:.2f} ({:d})\t{:.2f} ({:d})\t{:.2f} ({:d})\t{}\n'.format(100*sub_nan.shape[0]/Nnan, sub_nan.shape[0], 100*sub_big.shape[0]/big.shape[0], sub_big.shape[0], 100*sub_audit.shape[0]/xf.shape[0], sub_audit.shape[0], dac))

# looks like disproportionate amount of coriolis are nan or big - look into it
# sub = df[np.logical_or(df.diffGAIN.isnull(), df.diffGAIN >= 0.2)]
# sub = df[df.pyGAIN.isnull()]
sub = df[df.diffGAIN >= 1]
sub = sub[sub.DAC == 'coriolis']

sprof.set_dirs(argo_path='/Users/gordonc/Documents/data/Argo', woa_path='/Users/gordonc/Documents/data/WOA18')
syn = sprof(sub.WMO.iloc[0])
sf = syn.to_dataframe()
sf = sf[sf.CYCLE.isin([0,130])]

for flt in good.WMO.unique():
    syn = sprof(flt)
    # syn.clean()
    syn.calc_gains(ref='WOA')

    sub_df = good[good.WMO == flt]

    ff = syn.to_dataframe()
    ff = ff[ff.CYCLE.isin(sub_df.CYCLE)]
    ff = ff[ff.PRES < 30]

    fltWOA = syn.__WOAfloatref__
    WOA = syn.__WOAref__
    WOA = np.append(fltWOA, np.atleast_2d(WOA).T, axis=1)

    wf = pd.DataFrame(data=WOA, columns=['cycle', 'date', 'fltmean', 'fltstd', 'WOA'])
    wf = wf[wf.cycle.isin(sub_df.CYCLE)]

    sys.stdout.write('{}\n'.format(flt))
    sys.stdout.write('cycle\tpyG\tsageG\tflt\tWOA\n')
    for i,c in enumerate(sub_df.CYCLE):
        sys.stdout.write('{:d}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\n'.format(int(wf.cycle.iloc[i]), sub_df.pyGAIN.iloc[i], sub_df.sageGAIN.iloc[i], wf.fltmean.iloc[i], wf.WOA.iloc[i]))

    if plot_profiles:
        plt.plot(ff.O2Sat, ff.PRES)
        plt.show()