
from pathlib import Path
from .. import configure

class PathHandler:

    def __init__(self):

        config = configure.read_config()
        self.config = config
        self.ARGO_PATH = Path(config['argo_path']) if 'argo_path' in config.keys() else Path('./')
        self.NCEP_PATH = Path(config['ncep_path']) if 'ncep_path' in config.keys() else Path('./')
        self.WOA_PATH = Path(config['woa_path']) if 'woa_path' in config.keys() else Path('./')
    
    def set_dirs(self, **kwargs):
        '''
        Set local directories to look for Argo, WOA, and NCEP data.

        Args:
            argo_path (str or path-like): location of local Argo data
            ncep_data (str or path-like): location of local NCEP data
            woa_path (str or path-like): location of local World Ocean Atlas data
        '''

        self.ARGO_PATH = Path(kwargs.pop('argo_path')) if 'argo_path' in kwargs.keys() else self.ARGO_PATH
        self.NCEP_PATH = Path(kwargs.pop('ncep_path')) if 'ncep_path' in kwargs.keys() else self.NCEP_PATH
        self.WOA_PATH = Path(kwargs.pop('woa_path')) if 'woa_path' in kwargs.keys() else self.WOA_PATH

        if len(kwargs) > 0: # pragma: no cover
            raise ValueError(f'Invalid argument(s): {[k for k in kwargs.keys()]}')