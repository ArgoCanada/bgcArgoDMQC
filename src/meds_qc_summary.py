#!/usr/bin/python

from pathlib import Path
import glob

import numpy as np

from netCDF4 import Dataset

import matplotlib.pyplot as plt
import seaborn as sns
argo_path = Path('/Users/gordonc/Documents/data/Argo/meds/')
argo_glob = argo_path / '*'
float_paths = glob.glob(argo_glob.as_posix())

for fdir in float_paths[:1]:
    sprof = '{}_Sprof.nc'.format(fdir.split('\\')[-1])
    nc = Dataset(Path(fdir) / sprof, 'r')
    
    # oxygen and qc info
    DOXY = nc.variables['DOXY'][:]
    DOXY_ADJ = nc.variables['DOXY_ADJUSTED'][:]
    DOXY_QC = nc.variables['DOXY_QC'][:]
    DOXY_ADJ_QC = nc.variables['DOXY_ADJUSTED_QC'][:]
    
    # parameter info
    PARAMETER = nc.variables['PARAMETER']
    PARAMETER_DM = nc.variables['PARAMETER_DATA_MODE']

    # sci coeff info
    SCI_CAL_COEFF = nc.variables['SCIENTIFIC_CALIB_COEFFICIENT']
    SCI_CAL_COMMENT = nc.variables['SCIENTIFIC_CALIB_COMMENT']
    SCI_CAL_EQ = nc.variables['SCIENTIFIC_CALIB_EQUATION']
