#!/usr/bin/python

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

import sagepy

# ----------------------------------------------------------------------------
# example.py - script to demonstrate functionality of sagepy module
# ----------------------------------------------------------------------------
#
# Example script for demonstrating functionality that should run on most
# machines once provided with some path designations by the user. The script 
# is broken up into the following sections: 
#
# section 0: setup, designate paths
# section 1: download necessary data, if data already exists they won't be
# re-downloaded so requires no change form the user if paths are correct
# section 2: load and interpolate the data
# section 3: perform QC calculations with the data
# section 4: visualize the float and reference data, and computed QC info
#
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# section 0 - setting machine specific variables
# ----------------------------------------------------------------------------
# set the paths for data to be saved
# NOTE: sagepy uses python3's `pathlib` for handling paths, which means that
# all paths should be formatted in a unix-like way (omit C: on windows and
# use / instead of \) no matter what your operating system

# define where to save data - these directories must already exists or sagepy
# will throw errors - can make them by uncommenting lines 32-34
woa_path  = './WOA18'
ncep_path = './NCEP'
argo_path = './Argo'

woa_path  = '/Users/gordonc/Documents/data/WOA18'
ncep_path = '/Users/gordonc/Documents/data/NCEP'
argo_path = '/Users/gordonc/Documents/data/Argo/meds'

# uncomment these lines to make the above directories
# Path(woa_path).mkdir()
# Path(ncep_path).mkdir()
# Path(argo_path).mkdir()

# two argo floats - specific floats chosen to use as examples, but can be 
# altered with no performance change (as long as dac infomation corresponds
# to float numbers properly)
wmo_number = 4902481
dac_path = '/ifremer/argo/dac/meds'

# ----------------------------------------------------------------------------
# section 1 - downloading data
# ----------------------------------------------------------------------------
# download WOA and NCEP data necessary for calculating optode gain
sagepy.io.get_woa18('O2sat', local_path=woa_path) # oxygen % saturation
sagepy.io.get_ncep('pres', local_path=ncep_path) # surface air pressure in Pa

# download all data for argo floats
sagepy.io.get_argo(dac_path, [wmo_number], local_path=argo_path)

# ----------------------------------------------------------------------------
# section 2 - loading and interpolating reference data
# ----------------------------------------------------------------------------
# load in data from an argo float WITH NO in-air measurements
float_data = sagepy.argo(argo_path, wmo_number)

# make 'track' array with columns (time, lat, lon) to be used in interpolation
track = sagepy.track(float_data)

# get WOA data interpolated along track and associated depth levels
# NOTE: this function combines sagepy.io.load_woa_data() and 
# sagepy.interp.interp_woa_data() for convenience, but each can be called 
# individually if desired
z, woa_interp = sagepy.woa_to_float_track(track, 'O2sat', zlim=(0,1000), local_path=woa_path)
# put interpolated WOA data in dict for use later to calculate gains
ref_data = dict(z=z, WOA=woa_interp)

# get NCEP in-air data along track, again this is a convenience function and
# sagepy.io.load_ncep_data() and sagepy.interp.interp_ncep_data() can be called
# separately
ncep_pres_interp = sagepy.ncep_to_float_track('pres', track, local_path=ncep_path)/100

# ----------------------------------------------------------------------------
# section 3 - using the float and reference data to do QC
# ----------------------------------------------------------------------------

# calculate the gains based on WOA data, also returns the calculated surface
# data since WOA is depth resolved - specify limit here, default is 25dbar
woa_gains, surf_data = sagepy.calc_gain(float_data, ref_data, inair=False, zlim=25)

# gain calculation defaults to inair, don't need boolean
T = sagepy.get_var_by('TEMP_DOXY', 'TRAJ_CYCLE', float_data)
pH2O = sagepy.unit.pH2O(T)
ref_ppox = sagepy.unit.atmos_pO2(ncep_pres_interp[:-1], pH2O)
ncep_gains, inair_data = sagepy.calc_gain(float_data, ref_ppox)

# ----------------------------------------------------------------------------
# section 4 - plot the results
# ----------------------------------------------------------------------------
sdn = float_data['SDN']
fig, axes = plt.subplots(4,1,sharex=True)

g1 = sagepy.plt.float_woa_surface(sdn, surf_data[:,2], woa_interp[0,:], ax=axes[0])
g2 = sagepy.plt.gains(sdn, woa_gains, inair=False, ax=axes[1])

g3 = sagepy.plt.float_woa_surface(sdn[:-1], inair_data[:,2], ref_ppox, ax=axes[2])
g4 = sagepy.plt.gains(sdn[:-1], ncep_gains, inair=False, ax=axes[3])

plt.show()