#!/usr/bin/python

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns
# from pywaffle import Waffle

# summary comparison between bgcArgo and SAGE/DOXY audit
fn = Path('../data/argo_dmqc_comparison_20200818.h5')
df = pd.read_hdf(fn)
df['diffGAIN'] = df.pyGAIN - df.argoGAIN
df['absdiffGAIN'] = np.abs(df.pyGAIN - df.argoGAIN)

fig, axes = plt.subplots(1,2)
axes[0].plot(df.argoGAIN, df.pyGAIN, 'k.')
ll = np.min(axes[0].get_xlim() + axes[0].get_ylim())
ul = np.max(axes[0].get_xlim() + axes[0].get_ylim())
axes[0].plot((ll,ul), (ll,ul), 'k-')
axes[0].set_xlim((ll,ul))
axes[0].set_ylim((ll,ul))

axes[0].set_xlabel('$G_{Argo}$')
axes[0].set_ylabel('$G_{bgcArgo}$')

xf = df[df.diffGAIN < 1e9]
sns.distplot(xf.diffGAIN, kde=False, ax=axes[1])

axes[1].set_xlabel('$\Delta$G')

w, h = fig.get_figwidth(), fig.get_figheight()
fig.set_size_inches(w, h/2)
fig.tight_layout()
fig.savefig(Path('../figures/dmqc/scatter_and_dist_fulllims_20200818.png'), bbox_inches='tight', dpi=250)

fig, axes = plt.subplots(1,2)
axes[0].plot(df.argoGAIN, df.pyGAIN, 'k.')
ll = 0
ul = 10
axes[0].plot((ll,ul), (ll,ul), 'k-')
# axes[0].set_xlim((ll,ul))
# axes[0].set_ylim((ll,ul))

axes[0].set_xlabel('$G_{Argo}$')
axes[0].set_ylabel('$G_{bgcArgo}$')

xf = df[df.diffGAIN < 1e9]
sns.distplot(xf.diffGAIN, kde=False, ax=axes[1], bins=np.arange(0,2,0.05))

axes[1].set_xlabel('$\Delta$G')

w, h = fig.get_figwidth(), fig.get_figheight()
fig.set_size_inches(w, h/2)
fig.tight_layout()
fig.savefig(Path('../figures/dmqc/scatter_and_dist_medlims_20200818.png'), bbox_inches='tight', dpi=250)