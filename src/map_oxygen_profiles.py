#!/usr/bin/python

import numpy as np

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import seaborn as sns
import colorcet as cc

import bgcArgo as bgc

bp = bgc.get_index()
bp = bp[bp.parameters.notna()]
index = ['DOXY' in parameter_list for parameter_list in bp.parameters]
doxy = bp[index]

bin_size = 2
lat_bins = np.arange(-90, 9+bin_size, bin_size)
lon_bins = np.arange(-180, 18+bin_size, bin_size)

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
ax.set_global()
ax.coastlines()
ax.add_feature(cfeature.LAND)

ax.hist2d(doxy.longitude, doxy.latitude, bins=[lon_bins, lat_bins], cmap=cc.m_fire, transform=ccrs.PlateCarree())
plt.show()