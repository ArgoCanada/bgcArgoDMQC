
from pathlib import Path

def path(loc):

    local_absolute = Path(__file__).parent.absolute()
    loc = 'Argo/dac' if loc=='Argo' else loc

    print(local_absolute / loc)

    return local_absolute / loc