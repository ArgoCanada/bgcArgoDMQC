
from pathlib import Path

def path(loc):

    local_absolute = Path(__file__).parent.absolute()
    loc = 'Argo/dac' if loc=='Argo' else loc

    return local_absolute / loc