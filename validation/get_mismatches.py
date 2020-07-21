#!/usr/bin/python

import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns

# summary comparison between bgcArgo and SAGE/DOXY audit
fn = Path('../data/doxy_audit_vs_bgcArgo_py_comparison.csv')
df = pd.read_csv(fn)
df['diffGAIN'] = np.abs(df.pyGAIN - df.sageGAIN)

nan = df[df.diffGAIN.isnull()]
big = df[df.diffGAIN >= 0.2]

# by dac
for dac in nan.DAC.unique():
    sub_nan = nan[nan.DAC == dac]
    sub_big = big[big.DAC == dac]
    sys.stdout.write('{:d}\t{:d}\t{}\n'.format(sub_nan.shape[0], sub_big.shape[0], dac))

