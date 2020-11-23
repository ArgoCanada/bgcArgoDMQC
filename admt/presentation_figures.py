#!/usr/bin/python

from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import bgcArgoDMQC as bgc

sns.set(context='talk', style='ticks')

datapath = Path('/Users/gordonc/Documents/data/')
bgc.set_dirs(
    argo_path=datapath / 'Argo',
    ncep_path=datapath / 'NCEP',
    woa_path=datapath / 'WOA18'
)


### GAIN PLOT ###

# will load a synthetic (Sprof) file
syn  = bgc.sprof(4902481)
# calculate the gain using in-air data
gains = syn.calc_gains(ref='NCEP')
# plot the gains
g = syn.plot('gain', ref='NCEP')

# save figure
g.fig.savefig(Path('figures/gain_example.png'), bbox_inches='tight', dpi=350)
plt.close(g.fig)

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

syn.add_independent_data(label='CTD/Winkler', date=time, lat=lat, lon=lon, DOXY=doxy, TEMP=temp, PRES=pres)
fig, axes = syn.compare_independent_data()

fig.set_size_inches(8,6)
fig.savefig(Path('figures/independent_example.png'), bbox_inches='tight', dpi=350)
plt.close(fig)