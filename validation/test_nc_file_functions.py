#!/usr/bin/python

from pathlib import Path
from netCDF4 import Dataset

import numpy as np

import bgcArgoDMQC as bgc

myfile = Path('/Users/gordonc/Documents/data/Argo/meds/4902481/profiles/BR4902481_010.nc')
nc = Dataset(myfile, 'a')

exclude_dims = ['N_CALIB']
exclude_vars = [
    'SCIENTIFIC_CALIB_COMMENT', 
    'SCIENTIFIC_CALIB_EQUATION', 
    'SCIENTIFIC_CALIB_COEFFICIENT', 
    'SCIENTIFIC_CALIB_DATE',
    'PARAMETER',
    'DOXY_ADJUSTED',
    'DATA_MODE',
    'PARAMETER_DATA_MODE',
    'DOXY_ADJUSTED_QC',
    'PROFILE_DOXY_QC'
    ]
dims = [
    ('N_PROF', 'N_CALIB', 'N_PARAM', 'STRING256'),
    ('N_PROF', 'N_CALIB', 'N_PARAM', 'STRING256'),
    ('N_PROF', 'N_CALIB', 'N_PARAM', 'STRING256'),
    ('N_PROF', 'N_CALIB', 'N_PARAM', 'DATE_TIME'),
    ('N_PROF', 'N_CALIB', 'N_PARAM', 'STRING64'),
    ('N_PROF', 'N_LEVELS'),
    ('N_PROF'),
    ('N_PROF', 'N_PARAM'),
    ('N_PROF', 'N_LEVELS'),
    ('N_PROF')
]

sub_nc = bgc.io.copy_netcdf_except(
    myfile, 
    Path('/Users/gordonc/Documents/projects/bgcArgoDMQC/validation/BR4902481_010.nc'), 
    exclude_dims=exclude_dims, exclude_vars=exclude_vars
)

# create N_CALIB one greater than the previous one
sub_nc.createDimension('N_CALIB', size=nc.dimensions['N_CALIB'].size)
# create the variables that you removed, but with the new N_CALIB
for v,d in zip(exclude_vars, dims):
    sub_nc.createVariable(
        v, nc[v].datatype, fill_value=nc[v]._FillValue,
        dimensions=d
    )
    sub_nc[v].setncatts(nc[v].__dict__)

# make the first N_CALIB dimension the same as the last one
scientific_calib_comment = np.full(
    sub_nc['SCIENTIFIC_CALIB_COMMENT'].shape,
    sub_nc['SCIENTIFIC_CALIB_COMMENT']._FillValue,    
    dtype=sub_nc['SCIENTIFIC_CALIB_COMMENT'].datatype
    )

comment = 'A test comment written by Chris Gordon'
N = len(comment)
M = 256 - N
comment = comment + M*' '
comment = np.array(['{}'.format(let).encode('utf-8') for let in comment])
scientific_calib_comment[0, 0, 3, :] = comment
sub_nc['SCIENTIFIC_CALIB_COMMENT'][:] = scientific_calib_comment

scientific_calib_equation = np.full(
    sub_nc['SCIENTIFIC_CALIB_EQUATION'].shape,
    sub_nc['SCIENTIFIC_CALIB_EQUATION']._FillValue,    
    dtype=sub_nc['SCIENTIFIC_CALIB_EQUATION'].datatype
    )

equation = 'DOXY_ADJUSTED = G*DOXY - a test equation by Chris Gordon'
N = len(equation)
M = 256 - N
equation = equation + M*' '
equation = np.array(['{}'.format(let).encode('utf-8') for let in equation])
scientific_calib_equation[0, 0, 3, :] = equation
sub_nc['SCIENTIFIC_CALIB_EQUATION'][:] = scientific_calib_equation

scientific_calib_coefficient = np.full(
    sub_nc['SCIENTIFIC_CALIB_COEFFICIENT'].shape,
    sub_nc['SCIENTIFIC_CALIB_COEFFICIENT']._FillValue,    
    dtype=sub_nc['SCIENTIFIC_CALIB_COEFFICIENT'].datatype
    )

coeff = 'G = 1.04 - a made up coefficient for testing by Chris Gordon'
N = len(coeff)
M = 256 - N
coeff = coeff + M*' '
coeff = np.array(['{}'.format(let).encode('utf-8') for let in coeff])
scientific_calib_coefficient[0, 0, 3, :] = coeff
sub_nc['SCIENTIFIC_CALIB_COEFFICIENT'][:] = scientific_calib_coefficient

G = 1.04
doxy_adjusted = G*nc['DOXY'][:].data
sub_nc['DOXY_ADJUSTED'][:] = doxy_adjusted
doxy_adjusted_qc = nc['DOXY_QC'][:].data
doxy_adjusted_qc[doxy_adjusted_qc == b'3'] = b'2'
sub_nc['DOXY_ADJUSTED_QC'][:] = doxy_adjusted_qc
sub_nc['PARAMETER'][:] = nc['PARAMETER'][:].data
sub_nc['DATA_MODE'][:] = np.array([b'A'])
sub_nc['PROFILE_DOXY_QC'][:] = np.array([b'A'])
sub_nc['PARAMETER_DATA_MODE'][:] = np.array(['R','R','R','A','R'])

sub_nc.close()