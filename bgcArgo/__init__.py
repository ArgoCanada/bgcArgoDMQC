"""
package for quality control of BGC Argo floats
"""

from __future__ import absolute_import
from .core import *
from . import fplt
from . import unit
from . import io
from . import interp
from . import diagnostic

__all__ = ['fplt', 'unit', 'io', 'interp', 'diagnostic']

__author__ = ['Christopher Gordon <chris.gordon@dfo-mpo.gc.ca>']

__version__ = '0.2'
