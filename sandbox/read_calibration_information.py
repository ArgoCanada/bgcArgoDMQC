#!/usr/bin/python

from pathlib import Path
from netCDF4 import Dataset

import numpy as np

fn = Path('/Users/gordonc/Documents/data/Argo/aoml/3900345/profiles/BD3900345_004.nc')
nc = Dataset(fn)

eq    = nc.variables['SCIENTIFIC_CALIB_EQUATION']
coeff = nc.variables['SCIENTIFIC_CALIB_COEFFICIENT']
comm  = nc.variables['SCIENTIFIC_CALIB_COMMENT']

eqs = []
for row in np.squeeze(eq[:]):
    rval = ''
    for let in row.compressed():
        rval = rval + let.decode('UTF-8')
    eqs.append(rval)
