#!/usr/bin/python

import os, sys

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

try:
    import seaborn as sns
    sns.set(style='ticks', context='paper', palette='colorblind')
except:
    sys.stdout.write('Although the seaborn package is not required, it is recommended\n')

class pltClass:
    def __init__(self):
        self.info = 'Python qc package plt class'

def float_ncep_inair(sdn, flt, ncep, ax=None, legend=True):

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    ax.plot(sdn, flt, linewidth=2, label='Float')
    ax.plot(sdn, ncep, linewidth=2, label='NCEP')

    if legend:
        ax.legend(loc=3)

    mhr  = mdates.MonthLocator(interval=4)
    mihr = mdates.MonthLocator()
    fmt  = mdates.DateFormatter('%b %Y')

    ax.xaxis.set_major_locator(mhr)
    ax.xaxis.set_major_formatter(fmt)
    ax.xaxis.set_minor_locator(mihr)

    ax.set_ylabel('pO$_2$ (mbar)')

    for tick in ax.get_xticklabels():
        tick.set_rotation(45)

    g = pltClass()
    g.fig = fig
    g.ax  = ax

    return g

def float_woa_surface(sdn, flt, woa, ax=None, legend=True):

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    ax.plot(sdn, flt, linewidth=2, label='Float')
    ax.plot(sdn, woa, linewidth=2, label='WOA18')

    if legend:
        ax.legend(loc=3)

    mhr  = mdates.MonthLocator(interval=4)
    mihr = mdates.MonthLocator()
    fmt  = mdates.DateFormatter('%b %Y')

    ax.xaxis.set_major_locator(mhr)
    ax.xaxis.set_major_formatter(fmt)
    ax.xaxis.set_minor_locator(mihr)

    ax.set_ylabel('O$_2$ Saturation %')

    for tick in ax.get_xticklabels():
        tick.set_rotation(45)

    g = pltClass()
    g.fig = fig
    g.ax  = ax

    return g

def gains(sdn, gains, inair=True, ax=None, legend=True):

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    ax.plot(sdn, gains, 'o', markeredgewidth=0.5, markersize=5, markeredgecolor='grey', zorder=3, label='Gains')
    ax.axhline(np.nanmean(gains), color='k', linestyle='--', label='Mean = {:.2f}'.format(np.nanmean(gains)), zorder=2)
    ax.axhline(1.0, color='k', linestyle='-', linewidth=0.5, label=None,zorder=1)

    if legend:
        ax.legend(loc=3)

    mhr  = mdates.MonthLocator(interval=4)
    mihr = mdates.MonthLocator()
    fmt  = mdates.DateFormatter('%b %Y')

    ax.xaxis.set_major_locator(mhr)
    ax.xaxis.set_major_formatter(fmt)
    ax.xaxis.set_minor_locator(mihr)

    ax.set_ylabel('O$_2$ Gain (unitless)')

    for tick in ax.get_xticklabels():
        tick.set_rotation(45)

    g = pltClass()
    g.fig = fig
    g.ax  = ax

    return g

def gainplot(sdn, float_data, ref_data, gainvals, ref):

    fig, axes = plt.subplots(2,1,sharex=True)

    if ref == 'NCEP':

        g1 = float_ncep_inair(sdn, float_data, ref_data, ax=axes[0])
        g2 = gains(sdn, gainvals, inair=False, ax=axes[1])

    elif ref == 'WOA':

        g1 = float_woa_surface(sdn, float_data, ref_data, ax=axes[0])
        g2 = gains(sdn, gainvals, inair=False, ax=axes[1])

    return fig, axes
