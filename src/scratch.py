#!/usr/bin/python

from pathlib import Path

import numpy as np
from netCDF4 import Dataset

import bgcArgo as bgc

BRtraj = '/Users/gordonc/Documents/data/Argo/meds/4902481/4902481_BRtraj.nc'
nc = Dataset(Path(BRtraj), 'r')