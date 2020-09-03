# Argo Canada BGC Quality Control

## disclaimer

[![codecov](https://codecov.io/gh/ArgoCanada/BGC-QC/branch/master/graph/badge.svg)](https://codecov.io/gh/ArgoCanada/BGC-QC)

This code is in _very_ active development. Use of this code is available (encouraged even), but will likely throw errors, behave in undesired ways, etc. Submission of issues is also encouraged to help in development!

## installation

Available through the python package index <https://pypi.org/project/bgcArgo/>, install with:

```bash
pip install bgcArgo
```

Pull request has been submitted to get the recipe on conda-forge

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
from bgcArgo import profiles, sprof

# load data from profiles for two floats
flts = profiles([4902480, 6902905])
# calculate the dissolved oxygen gains
gains = profiles.calc_gain()
# visualize the oxygen gain QC step
fig, axes = profiles.plot('gains')

# load a synthetic profile
syn = sprof(4902480)
# plot a time vs. depth section for the top 500m
fig, ax = syn.plot('cscatter', 'DOXY', ylim=(0,500))
# plot the first 10 profiles for temperature, practical salinity,
# and adjusted oxygen
fig, axes = syn.plot('profiles', varlist=['TEMP','PSAL', 'DOXY_ADJUSTED'], Nprof=10)
```

## version history

0.1: April 20, 2020 - Initial creation

0.2: May 13, 2020 - Major change to how end user would use module, change to more object-oriented, create argo class

0.2.1: June 23, 2020 - pandas is now required, makes reading of global index significantly easier and more efficient

0.2.2: August 28, 2020 - remove pylab dependency (is part of matplotlib), built and uploaded to PyPI, build conda-forge recipe
