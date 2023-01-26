'''
# Argo Canada BGC Quality Control

## bgcArgoDMQC dependencies

- Must run on `python3.4` or higher, not supported on `python2.x` (uses [pathlib](https://docs.python.org/3/library/pathlib.html), introduced in python version 3.4)
- TEOS-10 package [gsw](https://teos-10.github.io/GSW-Python/), but will also work with the [seawater](https://pypi.org/project/seawater/) package, though it is deprecated in favor of gsw
- [netCDF4](https://pypi.org/project/netCDF4/) module for `.nc` files
- [pandas](https://pandas.pydata.org/) is required (and highly recommended for all your data science needs!)
- [seaborn](https://seaborn.pydata.org/) is recommended but not required, through there will be some reduced (non-essential) functionality
- [cmocean](https://matplotlib.org/cmocean/) is also recommended for nicer plots, but not required
'''

from __future__ import absolute_import

from . import configure
configure.check_config()

from . import plot
from . import unit
from . import util
from . import io
from . import interp

__all__ = ['plot', 'unit', 'util', 'io', 'interp', 'configure']

__author__ = ['Christopher Gordon <chris.gordon@dfo-mpo.gc.ca>']

__version__ = '0.2.13'

# check age of index file, or if it exists
if not io.index_exists():
    io.update_index()

from .core import *