#!/usr/bin/python

import os

import numpy as np
import pandas as pd

audit_file = '../data/external/DOXY_audit_022020.TXT'
gain_file  = '../data/external/DOXY_gains_022020.TXT'

af = pd.read_csv(audit_file, sep='\t', header=20)
# gf = pd.read_csv(gain_file, sep='\t', header=19)

# get only floats with MEDS as their DAC
meds = af.drop(af.index[af.DAC != 'meds'])
# how many different floats are in audit?
audit_floats = meds.WMO.unique()

DOWNLOAD=False
PLOT=False
SHALLOW=True

if DOWNLOAD:
    import ftplib
    # connect to dac ftp server
    url = 'ftp.ifremer.fr'
    dacdir = '/ifremer/argo/dac/meds'

    ftp = ftplib.FTP(url)
    ftp.login()
    ftp.cwd(dacdir)

    for flt in audit_floats:
        ftp.cwd(os.path.join(str(flt), 'profiles'))
        flt_df = meds.drop(meds.index[meds.WMO != flt])

        for cycle in flt_df.cycle:
            file = 'BR{}_{:03d}.nc'.format(flt,cycle)
            local_file = os.path.join('../data/audit/',str(flt),'profiles',file)
            lf = open(local_file, 'wb')
            ftp.retrbinary('RETR ' + file, lf.write)

        ftp.cwd('../../')

if PLOT:
    import matplotlib.pyplot as plt
    import seaborn as sns
    import cmocean.cm as cmo

    from netCDF4 import Dataset

    sns.set(style='ticks',palette='colorblind',context='paper')
    fig, axes = plt.subplots(1,len(audit_floats),sharey=True)
    zlim = 2000
    if SHALLOW:
        zlim = 150

    for i,flt in enumerate(audit_floats):
        flt_df = meds.drop(meds.index[meds.WMO != flt])
        zmin, zmax = flt_df.Z_raw.min(), flt_df.Z_raw.max()

        if zmin == zmax:
            zmin = zmin/2
            zmax = zmax + zmin
        for z,cycle in zip(flt_df.Z_raw,flt_df.cycle):
            file = 'BR{}_{:03d}.nc'.format(flt,cycle)
            local_file = os.path.join('../data/audit/',str(flt),'profiles',file)
            nc = Dataset(local_file,'r')

            P = nc.variables['PRES'][:]
            DOXY = nc.variables['DOXY'][:]

            P[P > zlim] = np.nan

            color = cmo.amp((z-zmin)/(zmax-zmin)+0.1)

            axes[i].plot(DOXY, P, 'o', markeredgecolor='k', markeredgewidth=0.05, color=color, markersize=5, alpha=0.5)

        axes[i].set_facecolor('#FAF0E6')
        axes[i].set_xlabel('Oxygen ($\mathregular{\mu}$mol kg$^{-1}$)')
        axes[i].set_title(str(flt),loc='left',weight='bold')
    axes[0].set_ylim((zlim,0))
    axes[0].set_ylabel('Depth (dbar)')

    w, h = fig.get_figwidth(), fig.get_figheight()
    fig.set_size_inches(w*1.25, h)

    filename = '../reports/figures/doxy_audit/profiles_zscore.png'
    if SHALLOW:
        filename = filename.replace('.png','_shallow.png')
    plt.savefig(filename, bbox_inches='tight', dpi=350)
    plt.close(fig)
