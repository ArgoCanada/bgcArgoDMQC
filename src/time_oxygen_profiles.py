#!/usr/bin/python

from pathlib import Path

import numpy as np
import pylab as pl

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns

import bgcArgo as bgc

sns.set(palette='colorblind', context='paper', style='ticks')

bp = bgc.get_index()
bp = bp[bp.parameters.notna()]
index = ['DOXY' in parameter_list for parameter_list in bp.parameters]
doxy = bp[index]

# drop na dates
nmissing = doxy.date.isna().sum()
doxy = doxy[doxy.date.notna()]

# interpret dates
bins = np.array([pl.datestr2num('{}{:02d}01'.format(y,m+1)) - 1 for y in range(2002, 2021) for m in range(12)])
bins = bins[:-5]
doxy['datenum'] = np.array([pl.datestr2num(str(d)[:-2]) for d in doxy.date])

ax = sns.distplot(doxy.datenum, bins=bins, kde=False, label='Monthly Profiles')
ax.set_xlabel('')
ax.set_title('Monthly Argo DOXY Profiles (N = {}, {} profiles missing date)'.format(doxy.shape[0]+nmissing, nmissing), loc='left')

ax2 = ax.twinx()
s, b = np.histogram(doxy.datenum, bins=bins)
s = np.cumsum(s)
bc = bins[:-1] + np.diff(bins)*0.5
ax2.plot(bc, s, lw=4, color=sns.color_palette('colorblind')[1])
ax.plot(np.nan, np.nan, lw=4, color=sns.color_palette('colorblind')[1], label='Cumulative Profiles')
ax.patch.set_visible(False)
ax.set_zorder(ax2.get_zorder() + 1)

ax.legend(loc=2)

for a in [ax, ax2]:
    a.set_xlim((bins[0], bins[-1]))

mihr = mdates.MonthLocator(range(1,13,3))
mhr  = mdates.YearLocator(2)
fmt  = mdates.DateFormatter('%Y')

ax.xaxis.set_major_locator(mhr)
ax.xaxis.set_major_formatter(fmt)
ax.xaxis.set_minor_locator(mihr)

ax.set_ylabel('Monthly Profiles')
ax2.set_ylabel('Cumulative Profiles')

plt.savefig(Path('../figures/argo_oxygen_profiles_times.png'), dpi=350, bbox_inches='tight')