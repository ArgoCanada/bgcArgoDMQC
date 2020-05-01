#!/usr/bin/python

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from scipy.io import loadmat

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
import seaborn as sns

import sagepy

sns.set(style='ticks', palette='colorblind', context='paper')

# -----------------------------------------------------------------------------
# Utility function definitions
# -----------------------------------------------------------------------------
def resid_plot(x, y1, y2, l1='$y_1$', l2='$y_2$', lr='$y_2 - y_1$', legend=True, xl='$x$', yl='$y$', ylr='residuals'):

    gs = GridSpec(2,3)
    fig = plt.figure()
    axes = [fig.add_subplot(gs[0,:2]), fig.add_subplot(gs[1,:2]),
            fig.add_subplot(gs[0,2]),  fig.add_subplot(gs[1,2])]

    axes[0].plot(x, y1, linewidth=2, label=l1)
    axes[0].plot(x, y2, linewidth=2, label=l2)
    axes[0].legend(loc=4)
    axes[0].set_ylabel(yl)

    axes[1].stem(x, y2-y1, label=lr, basefmt='k-')
    axes[1].legend(loc=4)
    axes[1].set_xlabel(xl)
    axes[1].set_ylabel(ylr)

    axes[2].plot(y1, y2, '.', zorder=2)

    li1 = np.min(np.append(np.abs(axes[2].get_xlim()), np.abs(axes[2].get_ylim())))
    li2 = np.max(np.append(np.abs(axes[2].get_xlim()), np.abs(axes[2].get_ylim())))
    l = [li1,li2]
    axes[2].plot(l,l,'k',zorder=1)
    axes[2].set_xlim(l)
    axes[2].set_ylim(l)

    axes[2].set_xlabel(l1)
    axes[2].set_ylabel(l2)

    sns.distplot(y2-y1, ax=axes[3], kde=False, label='$\mu = {:.3f}$'.format(np.nanmean(y2 - y1)))
    axes[3].set_xlabel(lr)
    axes[3].legend(loc=2, fontsize=8)

    return fig, axes

# -----------------------------------------------------------------------------
# Download MEDS oxygen floats
# -----------------------------------------------------------------------------
meds_float_list = 'meds_doxy_floats.csv'
df = pd.read_csv(meds_float_list)
df = df.drop(df.index[df.PROCESSED_BY == 'Coriolis'])

# remove leading Q from float numbers and any spaces
wmo_nums = np.array([n.strip('Q').strip(' ') for n in df.WMO_NUMBER])
ftp = sagepy.io.get_argo('/ifremer/argo/dac/meds', wmo_nums, local_path='/Users/gordonc/Documents/data/Argo/meds')

# -----------------------------------------------------------------------------
# Do gain calculation on an older float without in-air data
# -----------------------------------------------------------------------------
woa_path = '/Users/gordonc/Documents/data/WOA18'
local_path = '/Users/gordonc/Documents/data/Argo/meds'
wmo = '4900497'

data = sagepy.load_argo_data(local_path, wmo)
track = np.array([data['SDN'],data['LATITUDE'],data['LONGITUDE']]).T

xtrack, woa_track, woa_data = sagepy.load_woa_data(track, 'O2sat', zlim=(0,1000), local_path=woa_path)
woa_interp, wt, yrday = sagepy.interp_woa_data(xtrack, woa_track, woa_data)
z = woa_track[0]
# z, woa = sagepy.woa_to_float_track(track, 'O2sat', local_path=woa_path)

woa = dict(z=z, data=woa_interp)
gains, flt_surf, woa_surf = sagepy.calc_gain(data, woa, inair=False)

# -----------------------------------------------------------------------------
# Make analogous plots to SAGE-O2 GUI
# -----------------------------------------------------------------------------
sdn = data['SDN']
fig, axes = plt.subplots(2,1,sharex=True)

g1 = sagepy.plt.float_woa_surface(sdn, flt_surf[:,2], woa_surf, ax=axes[0])
g2 = sagepy.plt.gains(sdn, gains, inair=False, ax=axes[1])

fig.savefig(Path('./figures/sage_test/woa_float_gains.png'), dpi=350, bbox_inches='tight')
plt.close(fig)

# -----------------------------------------------------------------------------
# Compare python results to SAGE output
# -----------------------------------------------------------------------------

# mdict = loadmat('/Users/ChrisGordon/Documents/MATLAB/ARGO_PROCESSING/MFILES/GUIS/SAGE_O2Argo/test_data.mat')

# matlab_woa_data = mdict['d']
# matlab_woa_interp = mdict['d_interp']
# matlab_woa_surf = matlab_woa_interp[0,:]

# fig, axes = resid_plot(sdn, woa_surf, matlab_woa_surf, l1='WOA python', l2='WOA matlab SAGE', lr='matlab - python', xl='', yl='O2sat (%)')

mhr  = mdates.MonthLocator(interval=4)
mihr = mdates.MonthLocator()
fmt  = mdates.DateFormatter('%b %Y')

# axes[0].set_xticklabels([])
# axes[1].xaxis.set_major_locator(mhr)
# axes[1].xaxis.set_major_formatter(fmt)
# axes[1].xaxis.set_minor_locator(mihr)

# for tick in axes[1].get_xticklabels():
#     tick.set_rotation(45)

# plt.subplots_adjust(wspace=0.6, hspace=0.4)

# fig.savefig('../reports/figures/woa_py_matlab_compare.png', dpi=350, bbox_inches='tight')
# plt.close(fig)

mdict = loadmat(Path('/Users/gordonc/Documents/data/SAGE-O2/gui_data_4900297.mat'))
matlab_gains = np.squeeze(mdict['GAINS'])

fig, axes = resid_plot(sdn, gains, matlab_gains, l1='Gains python', l2='Gains matlab SAGE', lr='matlab - python', xl='', yl='Gain (unitless)')

axes[0].set_xticklabels([])
axes[1].xaxis.set_major_locator(mhr)
axes[1].xaxis.set_major_formatter(fmt)
axes[1].xaxis.set_minor_locator(mihr)

for tick in axes[1].get_xticklabels():
    tick.set_rotation(45)

plt.subplots_adjust(wspace=0.6, hspace=0.4)

fig.savefig('../reports/figures/gains_py_matlab_compare.png', dpi=350, bbox_inches='tight')
plt.close(fig)

matlab_surf_sat = np.squeeze(mdict['SURF_SAT'][:,1])

fig, axes = resid_plot(sdn, flt_surf[:,2], matlab_surf_sat, l1='Float python', l2='Float matlab SAGE', lr='matlab - python', xl='', yl='O2sat (%)')

axes[0].set_xticklabels([])
axes[1].xaxis.set_major_locator(mhr)
axes[1].xaxis.set_major_formatter(fmt)
axes[1].xaxis.set_minor_locator(mihr)

for tick in axes[1].get_xticklabels():
    tick.set_rotation(45)

plt.subplots_adjust(wspace=0.6, hspace=0.4)

fig.savefig('../reports/figures/float_py_matlab_compare.png', dpi=350, bbox_inches='tight')
plt.close(fig)
plt.show()