#!/usr/bin/python

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

# make new dataframe with categories depending on absolute deviation
counts = pd.DataFrame(dict(N=np.array([df.shape[0],
                            np.sum(np.isnan(df.pyGAIN)), 
                            df[np.logical_and(df.diffGAIN >= 0.01, df.diffGAIN < 0.05)].shape[0],
                            df[np.logical_and(df.diffGAIN >= 0.05, df.diffGAIN < 0.2)].shape[0],
                            df[df.diffGAIN >= 0.2].shape[0],
                            np.sum(np.logical_and(np.isinf(df.pyGAIN), np.isinf(df.sageGAIN))),
                            df[df.diffGAIN < 0.01].shape[0],]),
                        name=np.array(['Total', 'AD < 0.01', '0.01 <= AD < 0.05',
                                        '0.05 <= AD < 0.2', 'AD >= 0.2',
                                        'NaN valued', 'Both inf valued',])))
# calculate percent
counts[' '] = counts.N/counts.N.iloc[0]*100
# don't need total now
counts = counts[counts.name != 'Total']
colors = ['green', 'light green', 'yellow', 'red', 'grey', 'black']
palette = sns.xkcd_palette(colors)
sns.set(palette=palette, style='whitegrid', context='paper')

ax = counts.drop('N', axis=1).set_index('name').T.plot(kind='barh', stacked=True, width=1.0, legend=False)
ax.set_xlim((0,100))
ax.legend(loc=1, ncol=3, fontsize=8, bbox_to_anchor=(1.01, 1.6))
ax.set_aspect(10)
ax.set_title('$N={:d}$'.format(df.shape[0]), loc='left')
ax.xaxis.set_major_formatter(mtick.PercentFormatter())

plt.savefig('../figures/doxy_audit/DOXY_audit_comparison_breakdown.png', bbox_inches='tight', dpi=250)
plt.close()