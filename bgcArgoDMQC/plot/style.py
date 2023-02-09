
import matplotlib.pyplot as plt
cmocean_flag = True
try:
    import cmocean.cm as cmo
except:
    cmocean_flag = False

# default units for variables
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

# define colormaps
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
) if cmocean_flag else dict(
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