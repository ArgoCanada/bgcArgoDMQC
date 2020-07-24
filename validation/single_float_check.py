#!/usr/bin/python

import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import bgcArgo as bgc
from bgcArgo import sprof

flt_wmo = 5900245

fn = Path('../data/doxy_audit_vs_bgcArgo_py_comparison.csv')
df = pd.read_csv(fn)
df['diffGAIN'] = np.abs(df.pyGAIN - df.sageGAIN)

audit_file = Path('../data/DOXY_audit_070720.TXT')
xf = pd.read_csv(audit_file, sep='\t', header=25)
af = xf[xf.WMO == flt_wmo]

sprof.set_dirs(argo_path='/Users/gordonc/Documents/data/Argo', woa_path='/Users/gordonc/Documents/data/WOA18')
syn = sprof(flt_wmo)
syn.clean()
syn.calc_gains(ref='WOA')

ff = syn.to_dataframe()
ff = ff[ff.CYCLE.isin(af.cycle)]

gf = pd.DataFrame(dict(cycle = syn.CYCLE, gains = syn.gains))

xtrack, woa_track, woa_data = bgc.io.load_woa_data(syn.track, 'O2sat', zlim=(0,1000), local_path=syn.woa_path)