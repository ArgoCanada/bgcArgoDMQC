
import warnings

import pandas as pd

from .core import *
from .prof import prof
from .sprof import sprof
from .. import io

class float:
    '''
    Class that loads all relevant Argo data for a given float (wmo). 
    Loops bgcArgoDMQC.prof() for all available cycles, as well as loads
    synthetic and trajectory data through bgcArgoDMQC.sprof(). 
    '''

    def __init__(self, wmo, kind='bgc-b', direction=None):

        self.wmo = wmo
        
        # get local files
        glob = {
            'core':'*',
            'bgc-b':'B*',
            'bgc-s':'S*'
        }
        local_profile_files = get_local_profiles(io.Path.ARGO_PATH, wmo, glob=glob[kind])

        # get all available files from index
        index_profile_files = get_index_profiles(wmo, index=kind)

        # compare and warn if mismatch
        missing_files = list(set([Path(f).name for f in index_profile_files]) - set([Path(f).name for f in local_profile_files]))
        if len(missing_files) > 0:
            warning = f'There are {len(missing_files)} files in the index not found locally, consider running bgc.get_argo({wmo})'
            for f in missing_files:
                warning = warning + f'\n{f}'
            warning = warning + '\nOnly local files will be loaded'
            warnings.warn(warning)

        local_only_files = list(set([Path(f).name for f in local_profile_files]) - set([Path(f).name for f in index_profile_files]))
        if len(local_only_files) > 0:
            warning = f'There are {len(local_only_files)} files found only locally but not in the Argo index - has the data mode changed at the GDAC?'
            for f in local_only_files:
                warning = warning + f'\n{f}'
            warning = warning + '\nAll local files will be loaded'
            warnings.warn(warning)

        # load all prof() objects
        cycles = [f.as_posix().split('_')[-1].split('.')[0] for f in local_profile_files]
        direction = [cycle[-1] if cycle[-1] == 'D' else 'A' for cycle in cycles]
        cycles = [int(cycle[:-1]) if cycle[-1] == 'D' else int(cycle) for cycle in cycles]

        index = pd.MultiIndex.from_tuples([(a, b) for a, b in zip(cycles, direction)])
        profs = pd.Series([prof(file=fn) for fn in local_profile_files], index=index)

        self.profs = profs

        # load sprof() object
        self.sprof = sprof(wmo)

        # load all prof() data into a single dataframe, multi-indexed by cycle, nprof, nlevel