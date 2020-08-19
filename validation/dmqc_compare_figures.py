#!/usr/bin/python

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns
# from pywaffle import Waffle

# summary comparison between bgcArgo and SAGE/DOXY audit
fn = Path('../data/argo_dmqc_comparison_20200819.h5')
df = pd.read_hdf(fn)
df['diffGAIN'] = df.pyGAIN - df.argoGAIN
df['absdiffGAIN'] = np.abs(df.pyGAIN - df.argoGAIN)
df['diffMEANGAIN'] = df.pyMEANGAIN - df.argoGAIN
df['absdiffMEANGAIN'] = np.abs(df.pyMEANGAIN - df.argoGAIN)
df = df[df.absdiffMEANGAIN.notnull()]

fig, axes = plt.subplots(1,2)
axes[0].plot(df.argoGAIN, df.pyGAIN, 'k.')
ll = np.min(axes[0].get_xlim() + axes[0].get_ylim())
ul = np.max(axes[0].get_xlim() + axes[0].get_ylim())
axes[0].plot((ll,ul), (ll,ul), 'k-')
axes[0].set_xlim((ll,ul))
axes[0].set_ylim((ll,ul))

axes[0].set_xlabel('$G_{Argo}$')
axes[0].set_ylabel('$G_{bgcArgo}$')

sns.distplot(df.diffGAIN, kde=False, ax=axes[1])

axes[1].set_xlabel('$\Delta$G')

w, h = fig.get_figwidth(), fig.get_figheight()
fig.set_size_inches(w, h/2)
fig.tight_layout()
fig.savefig(Path('../figures/dmqc/scatter_and_dist_fulllims_20200819.png'), bbox_inches='tight', dpi=250)

fig, axes = plt.subplots(1,2)
axes[0].plot(df.argoGAIN, df.pyGAIN, 'k.')
ll = 0
ul = 2
axes[0].plot((ll,ul), (ll,ul), 'k-')
# axes[0].set_xlim((ll,ul))
# axes[0].set_ylim((ll,ul))

axes[0].set_xlabel('$G_{Argo}$')
axes[0].set_ylabel('$G_{bgcArgo}$')

sns.distplot(df.diffGAIN, kde=False, ax=axes[1], bins=np.arange(-2,2,0.05))

axes[1].set_xlabel('$\Delta$G')

w, h = fig.get_figwidth(), fig.get_figheight()
fig.set_size_inches(w, h/2)
fig.tight_layout()
fig.savefig(Path('../figures/dmqc/scatter_and_dist_medlims_20200819.png'), bbox_inches='tight', dpi=250)

wm = []
gv = []
gf = []
for w in df.WMO.unique():
    wm.append(w)
    gv.append(df[df.WMO == w].pyMEANGAIN.mean())
    gf.append(df[df.WMO == w].argoGAIN.mean())
wm = np.array(wm)
gv = np.array(gv)
gf = np.array(gf)

fig, ax = plt.subplots()
ax.plot(gf, gv, 'o', color=sns.color_palette('colorblind')[2])
ll = 0
ul = 2.0
ax.plot((ll,ul), (ll,ul), 'k-')
ax.set_xlim((ll,ul))
ax.set_ylim((ll,ul))

ax.set_xlabel('$G_{Argo}$')
ax.set_ylabel('$G_{bgcArgo}$')

fig.savefig(Path('../figures/dmqc/mean_scatter_20200819.png'), bbox_inches='tight', dpi=250)

fig, axes = plt.subplots(1,2)
axes[0].plot(df.argoGAIN, df.pyMEANGAIN, 'k.')
ll = 0
ul = 2
axes[0].plot((ll,ul), (ll,ul), 'k-')
axes[0].set_xlim((ll,ul))
axes[0].set_ylim((ll,ul))

axes[0].set_xlabel('$G_{Argo}$')
axes[0].set_ylabel('$G_{bgcArgo}$')

sns.distplot(df.diffMEANGAIN, kde=False, ax=axes[1], bins=np.arange(-2,2,0.05))

axes[1].set_xlabel('$\Delta$\mu_G')

w, h = fig.get_figwidth(), fig.get_figheight()
fig.set_size_inches(w, h/2)
fig.tight_layout()
fig.savefig(Path('../figures/dmqc/mean_scatter_and_dist_medlims_20200819.png'), bbox_inches='tight', dpi=250)

plt.close('all')