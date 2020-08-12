#!/usr/bin/python

from pathlib import Path
import pyowc as owc

from scipy.io import loadmat

# ----------------------------------------------------------------------------
# Example taken directly from pyowc github page, uses data that ships with
# the module so it should work no problem
# ----------------------------------------------------------------------------

# where the module data is saved on my machine
datapath = Path('/Users/gordonc/Documents/projects/external/argodmqc_owc/')
wmo = '3901960'
# creates an ordered dict that tells the module where to look for things
user_config = owc.configuration.load()
# make some changes to directories so it looks in the right places
for k in user_config: 
    if 'DIRECTORY' in k and 'data' in user_config[k][:5]:
        user_config[k] = str(datapath / user_config[k])

# print configuration information
print(owc.configuration.print_cfg(config))

# calculate mapped values needed for the analysis (verbatim from ipynb)
owc.calibration.update_salinity_mapping('', wmo, user_config)
# set the calseries parameteres for analysis and line fitting (verbatim from ipynb)
owc.configuration.set_calseries('', wmo, user_config)
# calculat the fit of each break and calibrate salinities (verbatim from ipynb)
owc.calibration.calc_piecewisefit('', wmo, user_config)

# produces a lot of plots
# owc.plot.dashboard('', wmo, user_config)

# ----------------------------------------------------------------------------
# Source data is all in .mat format, not traditional .nc format - is this
# because it is more raw (i.e. directly from float, I know APEX floats output
# both .txt and .mat files), or more processed (i.e. necessary data taken
# from the netcdf file) - investigate float file
# ----------------------------------------------------------------------------

# load the file
mdict = loadmat(datapath / 'data/float_source/3901960.mat')
# print out variable names
print(mdict.keys())

# ----------------------------------------------------------------------------
# Attempt at using the module to do something similar on actual Argo nc files
# ----------------------------------------------------------------------------
