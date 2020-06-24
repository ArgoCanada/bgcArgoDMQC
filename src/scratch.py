#!/usr/bin/python

from pathlib import Path

import numpy as np
from netCDF4 import Dataset

import bgcArgo as bgc
from argopy import DataFetcher as ArgoDataFetcher
argo_loader = ArgoDataFetcher()

# ds = argo_loader.region([-92.06, -62.06, 11.54, 41.54, 0, 1000]).to_xarray()
# ds = argo_loader.float(6902746).to_xarray().to_dataframe()

index = bgc.io.read_index()
wcs, matches = bgc.get_files(index, '/Users/gordonc/Documents/data/Argo/', [4902481, 4902416, 6902905])