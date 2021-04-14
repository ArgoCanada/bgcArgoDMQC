import numpy as np

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import seaborn as sns
sns.set(style='ticks', context='paper', palette='colorblind')

try:
    import cmocean.cm as cmo
    cmocean_flag = True
except:
    cmocean_flag = False

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

def var_cscatter(df, varname='DOXY', cmap=None, ax=None, ylim=(0,2000), clabel=None, vmin=None, vmax=None, **kwargs):
    # define colormaps
    if cmocean_flag:
        color_maps = dict(
            TEMP=cmo.thermal,
            TEMP_ADJUSTED=cmo.thermal,
            PSAL=cmo.haline, 
            PSAL_ADJUSTED=cmo.haline, 
            PDEN=cmo.dense,
            CHLA=cmo.algae,
            CHLA_ADJUSTED=cmo.algae,
            BBP700=cmo.matter,
            BBP700_ADJUSTED=cmo.matter,
            DOXY=cmo.ice,
            DOXY_ADJUSTED=cmo.ice,
            DOWNWELLING_IRRADIANCE=cmo.solar,
        )
    else:
        color_maps = dict(
            TEMP=plt.cm.inferno,
            TEMP_ADJUSTED=plt.cm.inferno,
            PSAL=plt.cm.viridis, 
            PSAL_ADJUSTED=plt.cm.viridis, 
            PDEN=plt.cm.cividis,
            CHLA=plt.cm.YlGn,
            CHLA_ADJUSTED=plt.cm.YlGn,
            BBP700=plt.cm.pink_r,
            BBP700_ADJUSTED=plt.cm.pink_r,
            DOXY=plt.cm.YlGnBu_r,
            DOXY_ADJUSTED=plt.cm.YlGnBu_r,
            DOWNWELLING_IRRADIANCE=plt.cm.magma,
        )

    if clabel is None:
        var_units = dict(
            TEMP='Temperature ({}C)'.format(chr(176)),
            TEMP_ADJUSTED='Temperature ({}C)'.format(chr(176)),
            PSAL='Practical Salinity', 
            PSAL_ADJUSTED='Practical Salinity', 
            PDEN='Potential Density (kg m${-3}$)',
            CHLA='Chlorophyll (mg m$^{-3}$',
            CHLA_ADJUSTED='Chlorophyll (mg m$^{-3}$',
            BBP700='$\mathsf{b_{bp}}$ (m$^{-1}$)',
            BBP700_ADJUSTED='$\mathsf{b_{bp}}$ (m$^{-1}$)',
            DOXY='Diss. Oxygen ($\mathregular{\mu}$mol kg$^{-1}$)',
            DOXY_ADJUSTED='Diss. Oxygen ($\mathregular{\mu}$mol kg$^{-1}$)',
            DOWNWELLING_IRRADIANCE='Downwelling Irradiance (W m$^{-2}$)',
        )
        clabel = var_units[varname]

    if cmap is None:
        cmap = color_maps[varname]
    
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    df = df.loc[df.PRES < ylim[1]+50]

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
        var_units = dict(
            TEMP='Temperature ({}C)'.format(chr(176)),
            TEMP_ADJUSTED='Temperature ({}C)'.format(chr(176)),
            PSAL='Practical Salinity', 
            PSAL_ADJUSTED='Practical Salinity', 
            PDEN='Potential Density (kg m$^{-3}$)',
            CHLA='Chlorophyll (mg m$^{-3}$',
            CHLA_ADJUSTED='Chlorophyll (mg m$^{-3}$',
            BBP700='$\mathsf{b_{bp}}$ (m$^{-1}$)',
            BBP700_ADJUSTED='$\mathsf{b_{bp}}$ (m$^{-1}$)',
            CDOM='CDOM (mg m$^{-3}$)',
            CDOM_ADJUSTED='CDOM (mg m$^{-3}$)',
            DOXY='Diss. Oxygen ($\mathregular{\mu}$mol kg$^{-1}$)',
            DOXY_ADJUSTED='Diss. Oxygen ($\mathregular{\mu}$mol kg$^{-1}$)',
            DOWNWELLING_IRRADIANCE='Downwelling Irradiance (W m$^{-2}$)',
        )
        xlabels = [var_units[v] for v in varlist]

    cm = plt.cm.gray

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

    for i,v in enumerate(varlist):
        for n in range(Nprof):
            subset_df = df.loc[df.CYCLE == CYCNUM[Ncycle-1 + n-1]]

            axes[i].plot(subset_df[v], subset_df[zvar], color=cm(CYCNUM[Ncycle-1 + n-1]/CYCNUM[-1]+0.05), **kwargs)
            
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
