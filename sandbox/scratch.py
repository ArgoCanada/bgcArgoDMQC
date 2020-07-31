#!/usr/bin/python

from pathlib import Path

import numpy as np
from netCDF4 import Dataset

import bgcArgo as bgc

# ds = argo_loader.region([-92.06, -62.06, 11.54, 41.54, 0, 1000]).to_xarray()
# ds = argo_loader.float(6902746).to_xarray().to_dataframe()

woa_path  = '/Users/gordonc/Documents/data/WOA18'
ncep_path = '/Users/gordonc/Documents/data/NCEP'
argo_path = '/Users/gordonc/Documents/data/Argo'
bgc.set_dirs(argo_path=argo_path, woa_path=woa_path, ncep_path=ncep_path)

flts = bgc.profiles([4902481, 6902905])
flts = flts.clean()
df = flts.to_dataframe()