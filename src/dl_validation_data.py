#!/usr/bin/python

from pathlib import Path
import numpy as np
import pandas as pd

import bgcArgo as bgc

def get_all_wmos():

    argo_path = Path('/Users/gordonc/Documents/data/Argo/')
    wmos = []
    dacs = list(argo_path.glob('*'))
    dacs.pop(-1)
    for dac in dacs:
        wmos = wmos + list(dac.glob('*'))

    wmos = [int(p.as_posix().split('/')[-1]) for p in wmos]
    
    return wmos

wmo = get_all_wmos()
dac = [bgc.get_dac(w) for w in wmo]

dacpath = '/ifremer/argo/dac'
fltpath = ['{}/{}/{}'.format(dacpath, d, w) for d, w in zip(dac, wmo)]

local_path = '/Users/gordonc/Documents/data/Argo'
# for imac
# local_path = '/Users/ChrisGordon/Desktop/argo/'
bgc.io.get_argo(fltpath, local_path=local_path, mission='B', mode='D')