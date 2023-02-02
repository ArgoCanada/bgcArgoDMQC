
from pathlib import Path

def path(loc):

    loc = 'Argo/dac' if loc=='Argo' else loc
    local_absolute = Path(__file__).parent.absolute() / 'data'

    return local_absolute / loc