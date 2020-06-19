# Argo Canada BGC Quality Control

`python` library of functions for quality controlling dissolved oxygen data.
Heavily based on the [SOCCOM BGC Argo QC methods](https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING) 
program in `matlab`, uses either 
[NCEP](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html) 
or [World Ocean Atlas](https://www.nodc.noaa.gov/OC5/woa18/) data to
calculate oxygen gains 
([*Johnson et al. 2015*](https://doi.org/10.1175/JTECH-D-15-0101.1)).

## bgcArgo dependencies

- Must run on `python3`, not supported on `python2.x` (uses [pathlib](https://docs.python.org/3/library/pathlib.html), introduced in python version 3.4)
- The [seawater](https://pypi.org/project/seawater/) package
- [netCDF4](https://pypi.org/project/netCDF4/) module for `.nc` files
- [pandas](https://pandas.pydata.org/) and [seaborn](https://seaborn.pydata.org/) are recommended but not required, through there will be some reduced (non-essential) functionality without pandas
- [cmocean](https://matplotlib.org/cmocean/) is also recommended for nicer plots, but not required

## version history

0.1: April 20, 2020 - Initial creation

0.2: May 13, 2020 - Major change to how end user would use module, change to more object-oriented, create argo class
