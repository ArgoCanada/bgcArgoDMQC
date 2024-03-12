
from .core import set_dirs
from .. import configure

# get a dict with with config info
config = configure.read_config()
# set the directories within the config file
dir_config = {k:v for k,v in config.items() if k in ['argo_path', 'woa_path', 'ncep_path']}
set_dirs(**dir_config)

from .core import *
from .oxygen import *
from .sprof import sprof
from .prof import prof