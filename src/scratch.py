#!/usr/bin/python

from pathlib import Path

import numpy as np
from netCDF4 import Dataset

import bgcArgo as bgc
from argopy import DataFetcher as ArgoDataFetcher
argo_loader = ArgoDataFetcher()

# ds = argo_loader.region([-92.06, -62.06, 11.54, 41.54, 0, 1000]).to_xarray()
ds = argo_loader.float(6902746).to_xarray().to_dataframe()