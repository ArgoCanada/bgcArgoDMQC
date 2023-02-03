
from pathlib import Path
from .. import configure

URL_DICT = {
    'ftp.ifremer.fr':'ftp.ifremer.fr', 
    'ifremer':'ftp.ifremer.fr', 
    'coriolis':'ftp.ifremer.fr', 
    'usgodae.org':'usgodae.org', 
    'godae':'usgodae.org', 
    'us':'usgodae.org'
}

URL_DIR_DICT = {
    'ftp.ifremer.fr':'/ifremer/argo/', 
    'usgodae.org':'/pub/outgoing/argo/', 
}

config = configure.read_config()
if 'default_url' in config.keys():
    url_name = config.pop('default_url')
    URL = URL_DICT[url_name]
else:
    URL = 'ftp.ifremer.fr'

def path(loc):

    loc = 'Argo/dac' if loc=='Argo' else loc
    local_absolute = Path(__file__).parent.absolute() / 'data'

    return local_absolute / loc