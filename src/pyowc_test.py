#!/usr/bin/python

from pathlib import Path
import pyowc as owc

from scipy.io import loadmat, savemat

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
user_config['FLOAT_PLOTS_FORMAT'] = 'pdf'

# print configuration information
print(owc.configuration.print_cfg(user_config))

# calculate mapped values needed for the analysis (verbatim from ipynb)
owc.calibration.update_salinity_mapping('', wmo, user_config)
# set the calseries parameteres for analysis and line fitting (verbatim from ipynb)
owc.configuration.set_calseries('', wmo, user_config)
# calculate the fit of each break and calibrate salinities (verbatim from ipynb)
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
print([k for k in mdict.keys()])

# ----------------------------------------------------------------------------
# Attempt at using the module to do something similar on actual Argo nc files
# am having trouble seeing if this module is even meant to handle .nc files, 
# so I am going to write a function to convert from .nc to .mat will the 
# valiables listed in the .mat file above. Will use sprof file. 
# ----------------------------------------------------------------------------

from netCDF4 import Dataset
from matplotlib.dates import datestr2num
import gsw

from bgcArgo import sprof
sprof.set_dirs(argo_path='/Users/gordonc/Documents/data/Argo')

def nc_to_owc_mat(ncfile):
    nc = Dataset(ncfile)
    export_dict = dict(
        PROFILE_NO = nc.variables['CYCLE_NUMBER'][:].data,
        LAT = nc.variables['LATITUDE'][:].data,
        LONG = nc.variables['LONGITUDE'][:].data,
        DATES = (nc.variables['JULD'][:].data + datestr2num('1950-01-01'))/365 + 1970, # note dates are in years which is gross
        PRES = nc.variables['PRES'][:].data,
        TEMP = nc.variables['TEMP'][:].data,
        PTMP = gsw.pt0_from_t(nc.variables['PSAL'][:].data, nc.variables['TEMP'][:].data, nc.variables['PRES'][:].data,),
        SAL = nc.variables['PSAL'][:].data,
    )

    return export_dict

wmo = 4902481 # a MEDS float
syn = sprof(wmo)
savemat(Path('../data/owc/{}.mat'.format(wmo)), nc_to_owc_mat(syn.__Sprof__))

# make changes to directories so it looks in the right place again
user_config['FLOAT_SOURCE_DIRECTORY'] = str(Path('../data/owc'))
user_config['FLOAT_CALIB_DIRECTORY'] = str(Path('../data/owc'))
user_config['FLOAT_MAPPED_DIRECTORY'] = str(Path('../data/owc'))
user_config['FLOAT_PLOTS_DIRECTORY'] = str(Path('../owc_plots'))

# try it
print(owc.configuration.print_cfg(user_config))

# calculate mapped values needed for the analysis (verbatim from ipynb)
owc.calibration.update_salinity_mapping('', str(wmo), user_config)
# set the calseries parameteres for analysis and line fitting (verbatim from ipynb)
owc.configuration.set_calseries('', str(wmo), user_config)
# calculate the fit of each break and calibrate salinities (verbatim from ipynb)
owc.calibration.calc_piecewisefit('', str(wmo), user_config)

owc.plot.dashboard('', str(wmo), user_config)