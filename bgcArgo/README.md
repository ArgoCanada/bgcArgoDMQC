# bgcArgo

`python` module created to emulate the [SOCCOM BGC Argo QC methods](https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING), as well as some additional functionality

## dependencies

- Must run on `python3`, not supported on `python2.x`
- The [seawater](https://pypi.org/project/seawater/) package
- [netCDF4](https://pypi.org/project/netCDF4/) module for `.nc` files
- [pandas](https://pandas.pydata.org/) and [seaborn](https://seaborn.pydata.org/) are recommended but not required

## version history

0.1: April 20, 2020 - Initial creation

0.2: May 13, 2020 - Major change to how end user would use module, change to more object-oriented, create argo class