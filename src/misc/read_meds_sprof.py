#!/usr/bin/python

import numpy as np

from netCDF4 import Dataset

filename = '../data/meds/4900494_Sprof.nc'
nc = Dataset(filename,'r')
