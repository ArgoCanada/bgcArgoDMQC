#!/usr/bin/python

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns
from pywaffle import Waffle

# summary comparison between bgcArgo and SAGE/DOXY audit
fn = Path('../data/doxy_audit_vs_bgcArgo_py_comparison_20200920.csv')
df = pd.read_csv(fn)
df['absdiffGAIN'] = np.abs(df.pyGAIN - df.sageGAIN)
df['diffGAIN'] = np.abs(df.pyGAIN - df.sageGAIN)

# make new dataframe with categories depending on absolute deviation
counts = pd.DataFrame(dict(N=np.array([df.shape[0],df[df.absdiffGAIN < 0.01].shape[0],
                            df[np.logical_and(df.absdiffGAIN >= 0.01, df.absdiffGAIN < 0.05)].shape[0],
                            df[np.logical_and(df.absdiffGAIN >= 0.05, df.absdiffGAIN < 0.2)].shape[0],
                            df[df.absdiffGAIN >= 0.2].shape[0],
                            np.sum(df.pyGAIN.isna()),
                            np.sum(np.logical_and(np.isinf(df.pyGAIN), np.isinf(df.sageGAIN)))]),
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

plt.savefig('../figures/doxy_audit/DOXY_audit_comparison_breakdown_20200919.png', bbox_inches='tight', dpi=250)
plt.close()

waf = plt.figure(FigureClass=Waffle, rows=5, values=counts[' '].values, labels=list(counts.name), colors=tuple(palette), 
    icons=['check-circle', 'check-circle', 'exclamation-triangle', 'times-circle', 'question-circle', 'infinity'], icon_legend=True, icon_size=10,
    title={'label': 'Absolute Deviation (AD) between python and SAGE-O2 gains\nfrom DOXY audit performed by Josh Plant, $N={:d}$'.format(df.shape[0]), 'loc': 'left'}, legend={'loc': 3, 'bbox_to_anchor': (-0.02, -0.5), 'ncol': 3, 'fontsize': 10})
# plt.savefig('../figures/doxy_audit/DOXY_audit_comparison_waffle_pct_20200730.png', bbox_inches='tight', dpi=250)
# plt.close()

waf = plt.figure(FigureClass=Waffle, rows=50, values=counts['N'].values, labels=list(counts.name), colors=tuple(palette), 
    icons=['check-circle', 'check-circle', 'exclamation-triangle', 'times-circle', 'question-circle', 'infinity'], icon_legend=True, icon_size=5,
    title={'label': 'Absolute Deviation (AD) between python and SAGE-O2 gains\nfrom DOXY audit performed by Josh Plant, $N={:d}$'.format(df.shape[0]), 'loc': 'left'}, 
    legend={'loc': 3, 'bbox_to_anchor': (-0.01, -0.135), 'ncol': 3, 'fontsize': 10})
plt.gcf().set_size_inches(10,10)
# plt.savefig('../figures/doxy_audit/DOXY_audit_comparison_waffle_20200919.png', bbox_inches='tight', dpi=250)
# plt.close()


fig, axes = plt.subplots(1,2)
axes[0].plot(df.sageGAIN, df.pyGAIN, 'k.')
ll = np.min(axes[0].get_xlim() + axes[0].get_ylim())
ul = np.max(axes[0].get_xlim() + axes[0].get_ylim())
axes[0].plot((ll,ul), (ll,ul), 'k-')
axes[0].set_xlim((ll,ul))
axes[0].set_ylim((ll,ul))

axes[0].set_xlabel('$G_{SAGE}$')
axes[0].set_ylabel('$G_{bgcArgo}$')

xf = df[df.diffGAIN.notnull()]
xf = xf[~np.isinf(xf.diffGAIN)]

sns.distplot(xf.diffGAIN, kde=False, ax=axes[1])

axes[1].set_xlabel('$\Delta$G')

w, h = fig.get_figwidth(), fig.get_figheight()
fig.set_size_inches(w, h/2)
fig.tight_layout()
fig.savefig(Path('../figures/doxy_audit/scatter_and_dist_fulllims_20200919.png'), bbox_inches='tight', dpi=250)

fig, axes = plt.subplots(1,2)
axes[0].plot(df.sageGAIN, df.pyGAIN, 'k.')
ll = 0
ul = 20
axes[0].plot((ll,ul), (ll,ul), 'k-')
axes[0].set_xlim((ll,ul))
axes[0].set_ylim((ll,ul))

axes[0].set_xlabel('$G_{SAGE}$')
axes[0].set_ylabel('$G_{bgcArgo}$')

sns.distplot(xf.diffGAIN, kde=False, ax=axes[1], bins=np.arange(-0.505, 0.505, 0.01))

axes[1].set_xlabel('$\Delta$G')

fig.set_size_inches(w, h/2)
fig.tight_layout()
fig.savefig(Path('../figures/doxy_audit/scatter_and_dist_medlims_20200919.png'), bbox_inches='tight', dpi=250)

fig, axes = plt.subplots(1,2)
axes[0].plot(df.sageGAIN, df.pyGAIN, 'k.')
ll = 0.5
ul = 2.5
axes[0].plot((ll,ul), (ll,ul), 'k-')
axes[0].set_xlim((ll,ul))
axes[0].set_ylim((ll,ul))

axes[0].set_xlabel('$G_{SAGE}$')
axes[0].set_ylabel('$G_{bgcArgo}$')

sns.distplot(xf.diffGAIN, kde=False, ax=axes[1], bins=np.arange(-0.205, 0.205, 0.01))

axes[1].set_xlabel('$\Delta$G')

fig.set_size_inches(w, h/2)
fig.tight_layout()
fig.savefig(Path('../figures/doxy_audit/scatter_and_dist_smalllims_20200919.png'), bbox_inches='tight', dpi=250)

plt.close('all')