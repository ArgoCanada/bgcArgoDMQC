#!/usr/bin/python

from pathlib import Path
from netCDF4 import Dataset

import bgcArgoDMQC as bgc

myfile = Path('/Users/gordonc/Documents/data/Argo/meds/4902481/4902481_Sprof.nc')
nc = Dataset(myfile, 'r')

exclude_dims = ['N_CALIB']
exclude_vars = [
    'SCIENTIFIC_CALIB_COMMENT', 
    'SCIENTIFIC_CALIB_EQUATION', 
    'SCIENTIFIC_CALIB_COEFFICIENT', 
    'PARAMETER'
    ]

sub_nc = bgc.io.copy_netcdf_except(myfile, 'tmp.nc', exclude_dims=exclude_dims, exclude_vars=exclude_vars)
