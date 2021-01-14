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
    'SCIENTIFIC_CALIB_DATE',
    'PARAMETER'
    ]

sub_nc = bgc.io.copy_netcdf_except(myfile, 'test_outfile.nc', exclude_dims=exclude_dims, exclude_vars=exclude_vars)

# create N_CALIB one greater than the previous one
sub_nc.createDimension('N_CALIB', size=nc.dimensions['N_CALIB'].size+1)
# create the variables that you removed, but with the new N_CALIB
for v in exclude_vars:
    sub_nc.createVariable(
        v, nc[v].datatype, fill_value=nc[v]._FillValue,
        dimensions=('N_PROF', 'N_CALIB', 'N_PARAM', 'STRING256')
    )