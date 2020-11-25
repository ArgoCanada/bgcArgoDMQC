#!/usr/bin/python

from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

from netCDF4 import Dataset

import bgcArgoDMQC as bgc

sns.set(context='talk', style='ticks')

datapath = Path('/Users/gordonc/Documents/data/')
bgc.set_dirs(
    argo_path=datapath / 'Argo',
    ncep_path=datapath / 'NCEP',
    woa_path=datapath / 'WOA18'
)


### GAIN PLOT ###
'''
# will load a synthetic (Sprof) file
syn  = bgc.sprof(4902481)
# calculate the gain using in-air data
gains = syn.calc_gains(ref='NCEP')
# plot the gains
g = syn.plot('gain', ref='NCEP')

# save figure
g.fig.savefig(Path('figures/gain_example.png'), bbox_inches='tight', dpi=350)
plt.close(g.fig)
'''
### INDEPENDENT DATA ###

files = list(Path('/Users/gordonc/Documents/data/Ship/PacOxyFloats/').glob('*chem*'))
fn = files[2]
df = pd.read_csv(fn, header=7)
wmo = int(fn.as_posix().split('_')[-1].split('.')[0])
syn = bgc.sprof(wmo)

with open(fn) as f:
    time = f.readline()[10:].rstrip()
    f.readline()
    lon = float(f.readline()[11:].rstrip())
    lat = float(f.readline()[10:].rstrip())

time = '17-Aug-2020 14:30:18'

mll  = df['Oxygen:Dissolved (mL/L)']
temp = df['Temperature:Secondary (\'deg)']
psal = df['Salinity:Bottle (PSS-78)']
pres = df['Pressure (decibar)']

doxy = bgc.unit.mL_per_L_to_umol_per_L(mll, temp)

syn.add_independent_data(label='Winkler', date=time, lat=lat, lon=lon, DOXY=doxy, PRES=pres)
syn.add_independent_data(label='CTD', date=time, lat=lat-1, lon=lon-1, TEMP=temp, PRES=pres)
fig, axes = syn.compare_independent_data()

for ax in axes[1:]:
    ax.set_title('')
    ax.set_ylabel('')
    ax.set_yticklabels([])

axes[1].set_title(f'Float #{wmo}')

fig.set_size_inches(10,6)
fig.savefig(Path('figures/independent_example.png'), bbox_inches='tight', dpi=350)

plt.close()

### RESPONSE TIME ###
'''
# load some time-resolved APEX data from GoM
fn = Path('/Users/gordonc/Documents/data/GoMRI/Sprof/f7939_Sprof.nc')
nc = Dataset(fn)

# extract one profile
ix = 3
time = nc.variables['MTIME'][:][ix,:].compressed()
pres = nc.variables['PRES'][:][ix,:].compressed()
doxy = nc.variables['DOXY'][:][ix,:].compressed()
temp = nc.variables['TEMP'][:][ix,:].compressed()

# do the time response correction
cd1 = bgc.correct_response_time_Tconst(time, doxy, 70)
cd2 = bgc.correct_response_time(time, doxy, temp, 125)

temp[pres > 300] = np.nan

fig, axes = plt.subplots(1, 2, sharey=True)

axes[0].plot(doxy, pres, label='Uncorrected')
axes[0].plot(cd1, pres, label='T-constant Corr.')
axes[0].plot(cd2, pres, label='T-dependent Corr.')
axes[1].plot(temp, pres)
axes[1].set_ylim((250,0))

axes[0].legend(loc=4, fontsize=10)

axes[0].set_xlabel('Diss. Oxygen ($\mathregular{\mu}$mol kg$^{-1}$)')
axes[1].set_xlabel('Temperature ({}C)'.format(chr(176)))
axes[0].set_ylabel('Pressure (dbar)')

fig.savefig(Path('figures/response_time_corr.png'), bbox_inches='tight', dpi=350)'''