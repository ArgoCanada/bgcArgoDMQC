# Argo Canada BGC Quality Control

## disclaimer

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/bgcargo/badges/installer/conda.svg)](https://conda.anaconda.org/conda-forge) [![Build Status](https://travis-ci.com/ArgoCanada/BGC-QC.svg?branch=master)](https://travis-ci.com/ArgoCanada/BGC-QC) [![codecov](https://codecov.io/gh/ArgoCanada/BGC-QC/branch/master/graph/badge.svg)](https://codecov.io/gh/ArgoCanada/BGC-QC)

This code is in _very_ active development. Use of this code is available (encouraged even), but will likely throw errors, behave in undesired ways, etc. Submission of issues is also encouraged to help in development!

## installation

The recommended install is through the conda-forge channel, via the command:

```bash
conda install -c conda-forge bgcargo
```

The package is also available through the python package index <https://pypi.org/project/bgcArgo/>, install with:

```bash
pip install bgcArgo
```

## general description

A `python` library of functions for quality controlling dissolved oxygen data.
Heavily based on the [SOCCOM BGC Argo QC methods](https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING)
program in `matlab`, uses either
[NCEP](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html)
or [World Ocean Atlas](https://www.nodc.noaa.gov/OC5/woa18/) data to
calculate oxygen gains
([*Johnson et al. 2015*](https://doi.org/10.1175/JTECH-D-15-0101.1)).

## bgcArgo dependencies

- Must run on `python3.4` or higher, not supported on `python2.x` (uses [pathlib](https://docs.python.org/3/library/pathlib.html), introduced in python version 3.4)
- TEOS-10 package [gsw](https://teos-10.github.io/GSW-Python/), but will also work with the [seawater](https://pypi.org/project/seawater/) package, though it is deprecated in favor of gsw
- [netCDF4](https://pypi.org/project/netCDF4/) module for `.nc` files
- [pandas](https://pandas.pydata.org/) is required (and highly recommended for all your data science needs!)
- [seaborn](https://seaborn.pydata.org/) is recommended but not required, through there will be some reduced (non-essential) functionality
- [cmocean](https://matplotlib.org/cmocean/) is also recommended for nicer plots, but not required

## basic functionality

Although functions in the `bgcArgo` module may be of use in other situations, the majority of the functionality is lies within two classes, `profiles` for typical profile files and `sprof` for synthetic profiles.

```python
import bgcArgo as bgc

# setup for your system - these directories need to already exist!
argo_path = 'your/argo/data/path' # where to save Argo data
ncep_path = 'your/ncep/data/path' # where to save NCEP reanalysis data
woa_path  = 'your/woa18/data/path' # where to save WOA data

# download the data - this can take some time depending on connection
# Argo
wmos = [4902480, 6902905]
dacs = [bgc.get_dac(w) for w in wmos]
dacpath = '/ifremer/argo/dac'
fltpath = ['{}/{}/{}'.format(dacpath, d, w) for d, w in zip(dacs, wmos)]
# NCEP
bgc.io.get_ncep('pres', local_path=ncep_path)
bgc.io.get_ncep('land', local_path=ncep_path)
# WOA
bgc.io.get_woa18('O2sat', local_path=woa_path)

# tell the package where to look for data
bgc.set_dirs(argo_path=argo_path, ncep_path=ncep_path, woa_path=woa_path)
# load data from profiles for two floats
flts = bgc.profiles(wmos)
# calculate the dissolved oxygen gains
gains = flts.calc_gains()
print(gains)

# load a synthetic profile
syn = bgc.sprof(4902480)
# plot a time vs. depth section for the top 500m
g1 = syn.plot('cscatter', varname='DOXY', ylim=(0,500))
# plot the first 10 profiles for temperature, practical salinity,
# and adjusted oxygen
g2 = syn.plot('profiles', varlist=['TEMP','PSAL', 'DOXY'], Ncycle=10, Nprof=10, ylim=(0,500))
```

The above code would produce the following plots:

![colorscale plot](https://raw.githubusercontent.com/ArgoCanada/BGC-QC/master/figures/example_p1.png | width=200)

![profiles plot](https://raw.githubusercontent.com/ArgoCanada/BGC-QC/master/figures/example_p2.png | width=200)

## version history

0.1: April 20, 2020 - Initial creation

0.2: May 13, 2020 - Major change to how end user would use module, change to more object-oriented, create argo class

0.2.1: June 23, 2020 - pandas is now required, makes reading of global index significantly easier and more efficient

0.2.2: August 28, 2020 - remove pylab dependency (is part of matplotlib), built and uploaded to PyPI, build conda-forge recipe

0.2.3 - 0.2.6: September 3, 2020 - updates to pass all checks on conda-forge pull request, updated on PyPI as well 
