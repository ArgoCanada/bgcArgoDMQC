#!/usr/bin/python

from pathlib import Path

import numpy as np
from netCDF4 import Dataset

import sagepy

wmo = 4902481
brtraj = '/Users/gordonc/Documents/data/Argo/meds/{}/{}_BRtraj.nc'.format(wmo,wmo)
sprof = brtraj.replace('BRtraj', 'Sprof')
bnc = Dataset(Path(brtraj), 'r')
snc = Dataset(Path(sprof), 'r')

fn = '/Users/gordonc/Documents/data/NCEP/pres.sfc.gauss.2020.nc'
ncep = Dataset(Path(fn), 'r')
