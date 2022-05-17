# Argo Canada BGC Quality Control

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/bgcargo/badges/installer/conda.svg)](https://anaconda.org/conda-forge/bgcargodmqc)
[![Build Status](https://travis-ci.com/ArgoCanada/bgcArgoDMQC.svg?branch=master)](https://travis-ci.com/ArgoCanada/bgcArgoDMQC)
[![Documentation Status](https://readthedocs.org/projects/bgcargodmqc/badge/?version=latest)](https://bgcargodmqc.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/ArgoCanada/bgcArgoDMQC/branch/master/graph/badge.svg)](https://codecov.io/gh/ArgoCanada/bgcArgoDMQC)

## installation

The recommended install is through the conda-forge channel, via the command:

```bash
conda install -c conda-forge bgcArgoDMQC
```

The package is also available through the python package index <https://pypi.org/project/bgcArgoDMQC/>, install with:

```bash
pip install bgcArgoDMQC
```
## dependencies

- Must run on `python3.4` or higher, not supported on `python2.x` (uses [pathlib](https://docs.python.org/3/library/pathlib.html), introduced in python version 3.4)
- TEOS-10 package [gsw](https://teos-10.github.io/GSW-Python/)
- [netCDF4](https://pypi.org/project/netCDF4/) module for `.nc` files
- [pandas](https://pandas.pydata.org/) is required
- [seaborn](https://seaborn.pydata.org/)
- [cmocean](https://matplotlib.org/cmocean/) is recommended for nicer plots, but not required

## setup

## general description

A `python` library of functions for quality controlling dissolved oxygen data.
Heavily based on the [SOCCOM BGC Argo QC methods](https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING)
program in `matlab`, uses either
[NCEP](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html)
or [World Ocean Atlas](https://www.nodc.noaa.gov/OC5/woa18/) data to
calculate oxygen gains
([*Johnson et al. 2015*](https://doi.org/10.1175/JTECH-D-15-0101.1)).

## basic functionality

## version history

0.1: April 20, 2020 - Initial creation

0.2: May 13, 2020 - Major change to how end user would use module, change to more object-oriented, create argo class

0.2.1: June 23, 2020 - pandas is now required, makes reading of global index significantly easier and more efficient

0.2.2: August 28, 2020 - remove pylab dependency (is part of matplotlib), built and uploaded to PyPI, build conda-forge recipe

0.2.3 - 0.2.6: September 3, 2020 - updates to pass all checks on conda-forge pull request, updated on PyPI as well

0.2.7 - 0.2.8: September 29, 2020 - re-spun for PyPI and PR to conda-feedstock

0.2.9: November 9, 2020 - name change to bgcArgoDMQC, various other updates over the past month
