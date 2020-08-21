#!/usr/bin/python

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import seaborn as sns
sns.set(style='ticks', context='paper', palette='colorblind')

def resid_plot(x, y1, y2, l1='$y_1$', l2='$y_2$', lr='$y_2 - y_1$', legend=True, xl='$x$', yl='$y$', ylr='residuals', axes=None):

    if axes is None:
        gs = GridSpec(2,3)
        fig = plt.figure()
        axes = [fig.add_subplot(gs[0,:2]), fig.add_subplot(gs[1,:2]),
                fig.add_subplot(gs[0,2]),  fig.add_subplot(gs[1,2])]
    else:
        fig = axes[0].get_figure()

    p = axes[0].plot(x, y1, linewidth=2, label=l1)
    axes[0].plot(x, y2, linewidth=2, label=l2, linestyle='--', color=p[0].get_color())
    axes[0].legend(loc=4, fontsize=4)
    axes[0].set_ylabel(yl)

    axes[1].stem(x, y2-y1, label=lr, basefmt='k-')
    # axes[1].legend(loc=4)
    axes[1].set_xlabel(xl)
    axes[1].set_ylabel(ylr)

    axes[2].plot(y1, y2, '.', zorder=2, alpha=0.5)

    li1 = np.min(np.append(np.abs(axes[2].get_xlim()), np.abs(axes[2].get_ylim())))
    li2 = np.max(np.append(np.abs(axes[2].get_xlim()), np.abs(axes[2].get_ylim())))
    l = [li1,li2]
    axes[2].plot(l,l,'k',zorder=1)
    # axes[2].set_xlim(l)
    # axes[2].set_ylim(l)

    axes[2].set_xlabel('sage')
    axes[2].set_ylabel('python')

    sns.distplot(y2-y1, ax=axes[3], kde=False, label='$\mu = {:.3f}$'.format(np.nanmean(y2 - y1)))
    axes[3].set_xlabel(lr)
    axes[3].legend(loc=1, fontsize=4)

    return fig, axes

df = pd.read_hdf(Path('../data/sage_gui_comparison.h5'))

fig = plt.figure()
gs = GridSpec(2,3)
axes = [fig.add_subplot(gs[0,:2]), fig.add_subplot(gs[1,:2]),
        fig.add_subplot(gs[0,2]),  fig.add_subplot(gs[1,2])]

# df.pyGAIN[df.pyGAIN > 2.4] = np.nan

for w in df.WMO.unique():
    xf = df[df.WMO == w]
    fig, axes = resid_plot(xf.CYCLE, xf.GAINS, xf.pyAIRGAIN, axes=axes, l1='{} SAGE'.format(w), l2='{} py'.format(w), xl='cycle', yl='gain', ylr='sage - py')

fig.tight_layout()
axes[2].set_xlim((0.9,1.3))
axes[2].set_ylim((0.9,1.3))
fig.savefig(Path('../figures/sage_test/sage_gui_comparison_air_cleaned.png'), bbox_inches='tight', dpi=250)
plt.close(fig)