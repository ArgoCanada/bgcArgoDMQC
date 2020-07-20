#!/usr/bin/python

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

# summary comparison between bgcArgo and SAGE/DOXY audit
fn = Path('../data/doxy_audit_vs_bgcArgo_py_comparison.csv')
df = pd.read_csv(fn)
df['diffGAIN'] = np.abs(df.pyGAIN - df.sageGAIN)

# make new dataframe with categories depending on absolute deviation
counts = pd.DataFrame(dict(N=np.array([df.shape[0],
                            np.sum(np.isnan(df.pyGAIN)), 
                            np.sum(np.logical_and(np.isinf(df.pyGAIN), np.isinf(df.sageGAIN))),
                            df[df.diffGAIN < 0.01].shape[0],
                            df[np.logical_and(df.diffGAIN >= 0.01, df.diffGAIN < 0.05)].shape[0],
                            df[np.logical_and(df.diffGAIN >= 0.05, df.diffGAIN < 0.2)].shape[0],
                            df[np.logical_and(df.diffGAIN >= 0.2, df.diffGAIN < 0.5)].shape[0],
                            df[df.diffGAIN >= 0.5].shape[0]]),
                        name=np.array(['Total', 'NaN valued', 'Both inf valued', 'AD < 0.01', '0.01 <= AD < 0.05',
                                        '0.05 <= AD < 0.2', '0.2 <= AD < 0.5', 'AD >= 0.5'])))
# calculate percent
counts['pct'] = counts.N/counts.N.iloc[0]*100
# double check that all data are captured by the categories
print(counts.pct.sum()-100)