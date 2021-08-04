from pathlib import Path
import sys

from .core import sprof
from . import io
from . import interp
from . import unit
from . import util
from . import plot

def simple_test(argo_path=None, ncep_path=None, woa_path=None, data_path=''):

    if any([argo_path is None, ncep_path is None, woa_path is None]):
        sys.stdout.write('You will be prompted to specify where you would like to store Argo, NCEP, and\nWOA data - please write relative or absolute paths as unix-like paths (i.e.\nusing ""/"" not ""\\"")\n')

    if all([argo_path is None, ncep_path is None, woa_path is None]) and data_path == '':
        data_path = input('If you would like to store all data in a single location, enter it here, or\nleave blank to specify each path individually: ')
    
    if data_path.strip() == '':
        if argo_path is None:
            argo_path = Path(input('Where would you like to store Argo data?: '))
        if ncep_path is None:
            ncep_path = Path(input('Where would you like to store NCEP data?: '))
        if woa_path is None:
            woa_path = Path(input('Where would you like to store WOA data?: '))
    else:
        data_path = Path(data_path)
        argo_path = data_path / 'Argo'
        ncep_path = data_path / 'NCEP'
        woa_path  = data_path / 'WOA18'

    for p in [argo_path, ncep_path, woa_path]:
        if not p.exists():
            p.mkdir()
    
    wmo = 4902480
    io.get_argo(wmo, local_path=argo_path, ftype='summary', overwrite=True)
    io.get_ncep('pres', local_path=ncep_path, years=[2019, 2020], overwrite=True)
    io.get_ncep('land', local_path=ncep_path, years=[2019, 2020], overwrite=True)
    io.get_woa18('O2sat', local_path=woa_path, overwrite=True, nfiles=2)
