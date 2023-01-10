import numpy as np

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.offsetbox import AnchoredText

import seaborn as sns
sns.set(style='ticks', context='paper', palette='colorblind')

try:
    import cmocean.cm as cmo
    cmocean_flag = True
except:
    cmocean_flag = False

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    carto_flag = True
except:
    carto_flag = False

from . import style
from .. import util

class pltClass:
    def __init__(self):
        self.__info__ = 'Python qc package plt class'

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
    g.fig  = fig
    g.axes = [ax]

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
    g.fig  = fig
    g.axes = [ax]

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
    g.fig  = fig
    g.axes = [ax]

    return g
    
def gainplot(sdn, float_data, ref_data, gainvals, ref):

    fig, axes = plt.subplots(2,1,sharex=True)

    if ref == 'NCEP':

        g1 = float_ncep_inair(sdn, float_data, ref_data, ax=axes[0])
        g2 = gains(sdn, gainvals, inair=False, ax=axes[1])

    elif ref == 'WOA':

        g1 = float_woa_surface(sdn, float_data, ref_data, ax=axes[0])
        g2 = gains(sdn, gainvals, inair=False, ax=axes[1])

    g = pltClass()
    g.fig  = fig
    g.axes = axes

    return g

def variable_color_scatter(df, varname='DOXY', cmap=None, ax=None, ylim=(0,2000), clabel=None, vmin=None, vmax=None, **kwargs):

    if clabel is None:
        clabel = style.var_units[varname]

    if cmap is None:
        cmap = style.color_maps[varname]
    
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    df = df[df.PRES < ylim[1]+50]

    if vmin is None:
        vmin = 1.05*df[varname].min()

    if vmax is None:
        vmax = 0.95*df[varname].max()

    im = ax.scatter(df.SDN, df.PRES, c=df[varname], s=50, cmap=cmap, vmin=vmin, vmax=vmax, **kwargs)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(clabel)
    ax.set_ylim(ylim)
    ax.invert_yaxis()

    ax.set_ylabel('Depth (dbar)')

    w, h = fig.get_figwidth(), fig.get_figheight()
    fig.set_size_inches(w*2, h)

    mhr  = mdates.MonthLocator(interval=4)
    mihr = mdates.MonthLocator()
    fmt  = mdates.DateFormatter('%b %Y')

    ax.xaxis.set_major_locator(mhr)
    ax.xaxis.set_major_formatter(fmt)
    ax.xaxis.set_minor_locator(mihr)

    g = pltClass()
    g.fig  = fig
    g.axes = [ax]
    g.cb = cb

    return g

def profiles(df, varlist=['DOXY'], Ncycle=1, Nprof=np.inf, zvar='PRES', xlabels=None, ylabel=None, axes=None, ylim=None, **kwargs):

    if xlabels is None:
        xlabels = [style.var_units[v] if v in style.var_units.keys() else '' for v in varlist]

    cm = plt.cm.gray_r

    if axes is None:
        fig, axes = plt.subplots(1, len(varlist), sharey=True)
        if len(varlist) == 1:
            axes = [axes]
    elif len(varlist) > 1:
        fig = axes[0].get_figure()
    else:
        fig = axes.get_figure()
        axes = [axes]

    if ylim is None:
        if zvar == 'PRES':
            ylim=(0,2000)
            if ylabel is None:
                ylabel = 'Pressure (dbar)'
        elif zvar == 'PDEN':
            ylim = (df.PDEN.min(), df.PDEN.max())
            if ylabel is None:
                ylabel = 'Density (kg m$^{-3}$)'

    df.loc[df[zvar] > ylim[1]*1.1] = np.nan

    CYCNUM = df.CYCLE.unique()
    greyflag = False
    if not 'color' in kwargs.keys():
        greyflag = True
    else:
        c = kwargs.pop('color')

    if Nprof > CYCNUM.shape[0]:
        Nprof = CYCNUM.shape[0]

    for i,v in enumerate(varlist):
        for n in range(Nprof):
            subset_df = df.loc[df.CYCLE == CYCNUM[Ncycle-1 + n-1]]

            if greyflag:
                c = cm(0.75*(CYCNUM[Ncycle-1 + n-1]/CYCNUM[-1])+0.25)
            axes[i].plot(subset_df[v], subset_df[zvar], color=c, **kwargs)
            
        axes[i].set_ylim(ylim[::-1])
        axes[i].set_xlabel(xlabels[i])

    subset_df = df.loc[df.CYCLE == CYCNUM[Ncycle-1]]
    date = mdates.num2date(subset_df.SDN.iloc[0]).strftime('%d %b, %Y')

    axes[0].set_ylabel(ylabel)
    if Nprof != 1:
        axes[0].set_title('Cyc. {:d}-{:d}, {}'.format(int(CYCNUM[Ncycle-1]), int(CYCNUM[Ncycle-1+Nprof-1]), date))
    else:
        axes[0].set_title('Cyc. {:d}, {}'.format(int(CYCNUM[Ncycle-1]), date))

    w, h = fig.get_figwidth(), fig.get_figheight()
    fig.set_size_inches(w*len(varlist)/3, h)

    g = pltClass()
    g.fig  = fig
    g.axes = axes

    return g

def qc_profiles(df, varlist=['DOXY'], Ncycle=1, Nprof=np.inf, zvar='PRES', xlabels=None, ylabel=None, axes=None, ylim=None, **kwargs):

    if xlabels is None:
        xlabels = [style.var_units[v] for v in varlist]

    if axes is None:
        fig, axes = plt.subplots(1, len(varlist), sharey=True)
        if len(varlist) == 1:
            axes = [axes]
    elif len(varlist) > 1:
        fig = axes[0].get_figure()
    else:
        fig = axes.get_figure()
        axes = [axes]

    if ylim is None:
        if zvar == 'PRES':
            ylim=(0,2000)
            if ylabel is None:
                ylabel = 'Pressure (dbar)'
        elif zvar == 'PDEN':
            ylim = (df.PDEN.min(), df.PDEN.max())
            if ylabel is None:
                ylabel = 'Density (kg m$^{-3}$)'

    df.loc[df[zvar] > ylim[1]*1.1] = np.nan

    CYCNUM = df.CYCLE.unique()

    if Nprof > CYCNUM.shape[0]:
        Nprof = CYCNUM.shape[0]

    groups = {'Good':[1,2,5], 'Probably Bad':[3], 'Bad':[4], 'Interpolated':[8]}
    colors = {'Good':'green', 'Probably Bad':'yellow', 'Bad':'red', 'Interpolated':'blue'}

    for i,v in enumerate(varlist):
        vqc = v + '_QC'
        for n in range(Nprof):
            subset_df = df.loc[df.CYCLE == CYCNUM[Ncycle-1 + n-1]]
            for k,f in groups.items():
                flag_subset_df = subset_df[subset_df[vqc].isin(f)]
                axes[i].plot(flag_subset_df[v], flag_subset_df[zvar], 'o', markeredgewidth=0.1, markeredgecolor='k', markerfacecolor=colors[k], **kwargs)
                
        axes[i].set_ylim(ylim[::-1])
        axes[i].set_xlabel(xlabels[i])

    subset_df = df.loc[df.CYCLE == CYCNUM[Ncycle-1]]
    date = mdates.num2date(subset_df.SDN.iloc[0]).strftime('%d %b, %Y')

    axes[0].set_ylabel(ylabel)
    if Nprof != 1:
        axes[0].set_title('Cyc. {:d}-{:d}, {}'.format(int(CYCNUM[Ncycle-1]), int(CYCNUM[Ncycle-1+Nprof-1]), date))
    else:
        axes[0].set_title('Cyc. {:d}, {}'.format(int(CYCNUM[Ncycle-1]), date))

    w, h = fig.get_figwidth(), fig.get_figheight()
    fig.set_size_inches(w*len(varlist)/3, h)

    g = pltClass()
    g.fig  = fig
    g.axes = axes

    return g

def compare_independent_data(df, plot_dict, meta_dict, fmt='*'):
    var_keys = []
    for label in plot_dict.keys():
        var_keys = var_keys + list(plot_dict[label].keys())
    var_keys = list(set(var_keys))
    var_keys.remove('PRES')

    meta_keys = []
    for label in meta_dict.keys():
        meta_keys = meta_keys + list(meta_dict[label].keys())
    meta_keys = set(meta_keys)

    meta_data_string = ''
    cyc = 1
    for label in meta_dict.keys():
        if label == ' ':
            label = 'Observation'

        if meta_dict[label]['date'] is None:
            cyc = 1
            dstr = 'N/A'
        else:
            cyc = util.cycle_from_time(meta_dict[label]['date'], df.SDN, df.CYCLE)
            dstr = mdates.num2date(meta_dict[label]['date']).strftime('%b %d, %Y')

        meta_data_string = meta_data_string + '{} date: {}\n'.format(label, dstr)
    meta_data_string = meta_data_string + 'Argo profile #{} date: {}'.format(cyc, mdates.num2date(df.SDN[df.CYCLE == cyc][0]).strftime('%b %d, %Y'))

    map_num = 0
    if 'lat' in meta_keys and 'lon' in meta_keys and carto_flag:
        map_num = 1 # change to 1 later, just broken right now
    
    nvar = len(set(var_keys))
    fig = plt.figure()
    ax_list = [fig.add_subplot(1, nvar+map_num, n+1) for n in range(nvar)]
    axes_dict = {v:ax for v, ax in zip(var_keys, ax_list)}
    if map_num > 0:
        ax_list.append(fig.add_subplot(1, nvar+map_num, nvar+1, projection=ccrs.PlateCarree()))
    
    ccount = 0
    fcol = sns.color_palette('colorblind')[0]
    clist = sns.color_palette('colorblind')[1:]
    
    for label in plot_dict.keys():
        pres = plot_dict[label].pop('PRES')

        if meta_dict[label]['date'] is None:
            cyc = 1
        else:
            cyc = util.cycle_from_time(meta_dict[label]['date'], df.SDN, df.CYCLE)

        varlist = list(plot_dict[label].keys())

        for v in varlist:
            profiles(df, varlist=[v], axes=axes_dict[v], Ncycle=cyc)
            axes_dict[v].plot(plot_dict[label][v], pres, fmt, label=None, color=clist[ccount])
        
        ax_list[0].plot(np.nan, np.nan, fmt, color=clist[ccount], label=label)
        ccount += 1

    ccount = 0
    if map_num > 0:
        c = []
        mx = ax_list[-1]
        dist_str = ''
        for label in meta_dict.keys():
            if meta_dict[label]['date'] is None:
                cyc = 1
            else:
                cyc = util.cycle_from_time(meta_dict[label]['date'], df.SDN, df.CYCLE)
            c1 = (meta_dict[label]['lat'], meta_dict[label]['lon'])
            c2 = (np.nanmean(df.LATITUDE[df.CYCLE == cyc]), np.nanmean(df.LONGITUDE[df.CYCLE == cyc]))
            c.append(c1)
            c.append(c2)
            mx.plot(c1[1], c1[0], fmt, transform=ccrs.PlateCarree(), label=label, color=clist[ccount])
            dist_str = dist_str + '\n{:.1f}km ({}) '.format(util.haversine(c1,c2), label)
            ccount += 1
        
            mx.plot(c2[1], c2[0], 'o', color=fcol, label=None, transform=ccrs.PlateCarree())
        mx.plot(np.nan, np.nan, 'o', color=fcol, label='Float {}'.format(df.WMO.iloc[0]))

        c = np.array(c)
        minlon, maxlon = np.nanmin(c[:,1]), np.nanmax(c[:,1])
        minlat, maxlat = np.nanmin(c[:,0]), np.nanmax(c[:,0])

        extent = [minlon, maxlon, minlat, maxlat]
        for i in range(len(extent)):
            if extent[i] < 0 and i % 2 == 0:
                extent[i] = 1.1*extent[i]
            elif extent[i] < 0 and i % 2 != 0:
                extent[i] = 0.9*extent[i]
            elif extent[i] > 0 and i % 2 == 0:
                extent[i] = 0.9*extent[i]
            elif extent[i] > 0 and i % 2 != 0:
                extent[i] = 1.1*extent[i]
        

        extent[2] = extent[2] - 6
        extent[3] = extent[3] + 6

        mx.set_extent(extent, crs=ccrs.PlateCarree())
        mx.legend(loc=4, bbox_to_anchor=(1.05, 1.0), fontsize=8)
        mx.add_feature(cfeature.COASTLINE)
        mx.add_feature(cfeature.BORDERS)
        mx.add_feature(cfeature.OCEAN)
        mx.add_feature(cfeature.LAND)
        mx.add_feature(cfeature.RIVERS)

        anc = AnchoredText('Distance between obs and float: ' + dist_str,
            loc=1, frameon=True, prop=dict(size=8))
        mx.add_artist(anc)
    
    if map_num == 0 and len(ax_list) > 1:
        for ax in ax_list[1:]:
            ax.set_title('')
            ax.set_ylabel('')
            ax.set_yticklabels([])
    elif map_num == 1 and len(ax_list) > 2:
        for ax in ax_list[1:-1]:
            ax.set_title('')
            ax.set_ylabel('')
            ax.set_yticklabels([])

    ax_list[0].legend(loc=2, fontsize=10)

    anc = AnchoredText(meta_data_string,
            loc=4, frameon=True, prop=dict(size=8))
    ax_list[0].add_artist(anc)

    fig.set_size_inches(10,6)

