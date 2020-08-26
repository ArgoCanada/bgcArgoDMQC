#!/usr/bin/python

from pathlib import Path
import sys

import numpy as np
from bgcArgo import sprof, profiles

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(context='paper', style='ticks', palette='colorblind')

# set up paths
data_path = Path('/Users/gordonc/Documents/data')
argo_path = data_path / 'Argo'
woa_path  = data_path / 'WOA18'
ncep_path = data_path / 'NCEP'

sprof.set_dirs(argo_path=argo_path, woa_path=woa_path, ncep_path=ncep_path)
profiles.set_dirs(argo_path=argo_path, woa_path=woa_path, ncep_path=ncep_path)

wmo = 4902481
# make synthetic profile object
syn = sprof(wmo)
# make profiles object for first 5 cycles only
prof = profiles(wmo, cycles=list(range(1,6)))

# calculate gains using sprof() object via NCEP and WOA
syn_air_gains = syn.calc_gains(ref='NCEP')[syn.CYCLE < 6]
syn_sat_gains = syn.calc_gains(ref='WOA')[syn.CYCLE < 6]

# calculate gains using profiles() object
# this doesn't work - no in air not in traj files?
# prof_air_gains = prof.calc_gains(ref='NCEP')
prof_sat_gains = prof.calc_gains(ref='WOA')

print('Synthetic mean in-air gain (first 5 cycles): {:.2f}'.format(np.mean(syn_air_gains)))
print('Synthetic WOA %-saturation mean gain: {:.2f}'.format(np.mean(syn_sat_gains)))
print('Profiles WOA %-saturation mean gain: {:.2f}'.format(np.mean(prof_sat_gains)))

fig,ax = plt.subplots()
ax.plot(syn_sat_gains, prof_sat_gains, 'o', zorder=2, label='Syn. WOA vs. Prof WOA')
ax.plot(syn_sat_gains, syn_air_gains, 'o', zorder=3, label='Syn. WOA vs. Syn. in-air')
yl = ax.get_ylim()
ax.plot(yl, yl, 'k-', zorder=1)
ax.set_ylim(yl)
ax.set_xlim(yl)

ax.set_xlabel('Synthetic WOA gains')
ax.set_ylabel('Profiles WOA gains OR Synthetic in-air gains')

ax.legend(loc=4)

plt.show()