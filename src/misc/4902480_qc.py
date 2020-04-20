#!/usr/bin/python

import os

import numpy as np
import pandas as pd

from netCDF4 import Dataset

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import cmocean.cm as cmo

sns.set(style='ticks',palette='colorblind',context='paper')

audit_file = '../data/external/DOXY_audit_meds_041620.TXT'
meds = pd.read_csv(audit_file, sep='\t', header=20)

# pick our 4902480 to look at
flt = 4902480
df = meds.drop(meds.index[meds.WMO != flt])

SHALLOW = False
append = ''

if True:
    """plot all doxy profiles for float, and audit flagged profiles"""
    fig, axes = plt.subplots(1,2,sharey=True)
    zlim = 2200
    if SHALLOW:
        zlim = 200
        append = '_shallow'

    cmap = sns.dark_palette('light blue', input='xkcd', as_cmap=True) # define the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist, cmap.N)

    bounds = np.arange(1, 32, 1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    for i in range(31):
        file = 'BR{}_{:03d}.nc'.format(flt,i+1)
        local_file = os.path.join('../data/meds/',str(flt),'profiles',file)
        nc = Dataset(local_file,'r')

        P = np.squeeze(nc.variables['PRES'][:])
        DOXY = np.squeeze(nc.variables['DOXY'][:])

        if i == 0:
            print(DOXY)

        axes[0].plot(DOXY, P, linewidth=2, color=cmap(i/30))


    fake_cbp = axes[0].scatter([np.nan,np.nan], [np.nan,np.nan], c=[10,15], vmin=1, vmax=31, cmap=cmap, norm=norm)
    cb = plt.colorbar(fake_cbp,ax=axes[0], cmap=cmap, norm=norm,
        spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i', label='Cycle #')

    zmin, zmax = 5, 20
    for z,cycle in zip(df.Z_raw,df.cycle):
        file = 'BR{}_{:03d}.nc'.format(4902480,cycle)
        local_file = os.path.join('../data/meds/',str(flt),'profiles',file)
        nc = Dataset(local_file,'r')

        P = np.squeeze(nc.variables['PRES'][:])
        DOXY = np.squeeze(nc.variables['DOXY'][:])

        P[P > zlim] = np.nan

        color = cmo.amp((z-zmin)/(zmax-zmin)+0.1)

        axes[1].plot(DOXY, P, 'o', markeredgecolor='k', markeredgewidth=0.1, color=color, markersize=5, alpha=0.5)
    fake_cbp = axes[1].scatter([np.nan,np.nan], [np.nan,np.nan], c=[10,15], vmin=zmin, vmax=zmax, cmap=cmo.amp)
    cb = plt.colorbar(fake_cbp,ax=axes[1],label='Z-score')

    axes[0].set_xlabel('Oxygen ($\mathregular{\mu}$mol kg$^{-1}$)')
    axes[1].set_xlabel('Oxygen ($\mathregular{\mu}$mol kg$^{-1}$)')
    axes[0].set_title('{} All Profiles'.format(flt),loc='left',weight='bold')
    axes[1].set_title('{} Audit Profiles'.format(flt),loc='left',weight='bold')
    axes[0].set_ylim((zlim,0))
    axes[0].set_ylabel('Depth (dbar)')

    fig.savefig('../reports/figures/4902480_zscore_profiles{}.png'.format(append), dpi=350, bbox_inches='tight')

if True:
    """plot doxy, phase for flagged floats"""
    fig, axes = plt.subplots(1, 4, sharey=True)
    zlim = 2200
    if SHALLOW:
        zlim = 200
        app = '_shallow'

    cmap = sns.dark_palette('light blue', input='xkcd', as_cmap=True) # define the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist, cmap.N)

    bounds = np.arange(1, df.shape[0], 1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    lbls = ['C1 Phase', 'C2 Phase', 'DOXY', 'TEMP_DOXY']

    for i,cycle in enumerate(df.cycle):
        file = 'BR{}_{:03d}.nc'.format(4902480,cycle)
        local_file = os.path.join('../data/meds/',str(flt),'profiles',file)
        nc = Dataset(local_file,'r')

        PRES = np.squeeze(nc.variables['PRES'][:])
        PHASE1 = np.squeeze(nc.variables['C1PHASE_DOXY'])
        PHASE2 = np.squeeze(nc.variables['C2PHASE_DOXY'])
        DOXY = np.squeeze(nc.variables['DOXY'][:])
        # TEMP = np.squeeze(nc.variables['TEMP'][:])
        TEMP_DOXY = np.squeeze(nc.variables['TEMP_DOXY'][:])
        # PSAL = np.squeeze(nc.variables['PSAL'][:])

        vars = [PHASE1, PHASE2, DOXY, TEMP_DOXY,]


        for v,ax in zip(vars, axes):
            ax.plot(v, PRES, linewidth=2, color=cmap(i/bounds.shape[0]))

    for l,ax in zip(lbls,axes):
        ax.set_xlabel(l)

    axes[0].set_ylim((zlim,0))
    axes[0].set_ylabel('Depth (dbar)')

    w, h = fig.get_figwidth(), fig.get_figheight()
    fig.set_size_inches(w*1.33, h)

    fig.savefig('../reports/figures/4902480_doxy_profiles{}.png'.format(append), dpi=350, bbox_inches='tight')

if True:
    """plot physics for flagged floats"""
    fig, axes = plt.subplots(1, 2, sharey=True)
    zlim = 2200
    if SHALLOW:
        zlim = 200
        app = '_shallow'

    cmap = sns.dark_palette('light blue', input='xkcd', as_cmap=True) # define the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist, cmap.N)

    bounds = np.arange(1, df.shape[0], 1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    lbls = ['TEMP', 'PSAL']

    for i,cycle in enumerate(df.cycle):
        file = 'R{}_{:03d}.nc'.format(4902480,cycle)
        local_file = os.path.join('../data/meds/',str(flt),'profiles',file)
        nc = Dataset(local_file,'r')

        PRES = np.squeeze(nc.variables['PRES'][:])
        TEMP = np.squeeze(nc.variables['TEMP'][:])
        PSAL = np.squeeze(nc.variables['PSAL'][:])

        vars = [TEMP, PSAL]

        print(nc.variables['PSAL_QC'][:])

        for v,ax in zip(vars, axes):
            ax.plot(v, PRES, linewidth=2, color=cmap(i/bounds.shape[0]))

    for l,ax in zip(lbls,axes):
        ax.set_xlabel(l)

    axes[0].set_ylim((zlim,0))
    axes[0].set_ylabel('Depth (dbar)')

    w, h = fig.get_figwidth(), fig.get_figheight()
    fig.set_size_inches(w*0.66, h)

    fig.savefig('../reports/figures/4902480_physics_profiles{}.png'.format(append), dpi=350, bbox_inches='tight')

plt.close('all')
