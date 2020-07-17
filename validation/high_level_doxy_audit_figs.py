#!/usr/bin/python

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

fn = Path('../data/doxy_audit_vs_bgcArgo_py_comparison.csv')
df = pd.read_csv(fn)
df['diffGAIN'] = np.abs(df.pyGAIN - df.sageGAIN)

counts = dict(total=df.shape[0], 
     nan=np.sum(np.isnan(df.pyGAIN)), 
     within_one_onehundredth=df[df.diffGAIN < 0.01].shape[0])
counts = pd.DataFrame(dict(N=np.array([df.shape[0],
                            np.sum(np.isnan(df.pyGAIN)), 
                            df[df.diffGAIN < 0.01].shape[0],
                            df[np.logical_and(df.diffGAIN >= 0.01, df.diffGAIN < 0.2)].shape[0],
                            df[np.logical_and(df.diffGAIN >= 0.2, df.diffGAIN < 0.5)].shape[0],
                            df[df.diffGAIN >= 0.5].shape[0]]),
                        name=np.array(['Total', 'NaN valued', 'AD < 0.01', '0.01 <= AD < 0.2',
                                        '0.2 <= AD < 0.5', 'AD >= 0.5'])))