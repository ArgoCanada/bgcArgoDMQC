'''
# Argo Canada BGC Quality Control

## bgcArgo dependencies

- Must run on `python3.4` or higher, not supported on `python2.x` (uses [pathlib](https://docs.python.org/3/library/pathlib.html), introduced in python version 3.4)
- TEOS-10 package [gsw](https://teos-10.github.io/GSW-Python/), but will also work with the [seawater](https://pypi.org/project/seawater/) package, though it is deprecated in favor of gsw
- [netCDF4](https://pypi.org/project/netCDF4/) module for `.nc` files
- [pandas](https://pandas.pydata.org/) is required (and highly recommended for all your data science needs!)
- [seaborn](https://seaborn.pydata.org/) is recommended but not required, through there will be some reduced (non-essential) functionality
- [cmocean](https://matplotlib.org/cmocean/) is also recommended for nicer plots, but not required

## version history

0.1: April 20, 2020 - Initial creation

0.2: May 13, 2020 - Major change to how end user would use module, change to more object-oriented, create argo class

0.2.1: June 23, 2020 - pandas is now required, makes reading of global index significantly easier and more efficient

0.2.2: August 28, 2020 - remove pylab dependency (is part of matplotlib), built and uploaded to PyPI, build conda-forge recipe

0.2.3 - 0.2.6: September 3, 2020 - updates to pass all checks on conda-forge pull request, updated on PyPI as well

0.2.7 - 0.2.8: September 29, 2020 - re-spun for PyPI and PR to conda-feedstock
'''

from __future__ import absolute_import
from .core import *
from . import fplt
from . import unit
from . import util
from . import io
from . import interp
from . import diagnostic
from . import calibration

__all__ = ['fplt', 'unit', 'util', 'io', 'interp', 'diagnostic', 'calibration']

__author__ = ['Christopher Gordon <chris.gordon@dfo-mpo.gc.ca>']

__version__ = '0.2.8'

# check age of index file, or if it exists
if not io.index_exists():
    io.update_index()
io.check_index()
