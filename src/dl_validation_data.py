#!/usr/bin/python

import numpy as np
import pandas as pd

import bgcArgo

wmo = [ 5904179, 5904682, 4902414, 4901779,
        4902481, 1901879, 6902905, 6903237, 6903273, 6902805 ]
dac = 2*['aoml'] + 3*['meds'] + ['bodc'] + 4*['coriolis']

wmo = [3900407]
dac = ['aoml']

dacpath = '/ifremer/argo/dac'
fltpath = ['{}/{}/{}'.format(dacpath, d, w) for d, w in zip(dac, wmo)]

local_path = '/Users/gordonc/Documents/data/Argo'
# for imac
# local_path = '/Users/ChrisGordon/Desktop/argo/'
bgcArgo.io.get_argo(fltpath, local_path=local_path)