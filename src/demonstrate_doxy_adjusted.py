#!/usr/bin/python

from bgcArgo import profiles

from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='ticks', palette='colorblind')

datapath = Path('/Users/gordonc/Documents/data/Argo')
profiles.set_dirs(argo_path=datapath)
flts = profiles([4902481, 6902905])
flts.clean()

df = flts.to_dataframe()

fig, ax = plt.subplots()
sns.scatterplot(x='DOXY', y='PRES', style='WMO', hue='CYCLE', data=df, ax=ax, linewidth=0.1, alpha=0.5)
ax.set_xlim((80,340))
ax.set_ylim((1000,0))
plt.show()
