"""
package for quality control of BGC Argo floats
"""

from __future__ import absolute_import
from pathlib import Path
from .core import *
from . import unit
from . import plt
from . import io

__all__ = ['plt', 'unit', 'io']

__author__ = ['Christopher Gordon <chris.gordon@dfo-mpo.gc.ca>']

__version__ = '0.2'

set_dirs()

if str(Path.home()) == 'C:\\Users\\gordonc':
    argo_path = '/Users/gordonc/Documents/data/Argo/meds'
    woa_path  = '/Users/gordonc/Documents/data/WOA18'
    ncep_path  = '/Users/gordonc/Documents/data/NCEP'

    set_dirs(argo_path=argo_path, woa_path=woa_path, ncep_path=ncep_path)