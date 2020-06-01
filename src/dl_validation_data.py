#!/usr/bin/python

import numpy as np
import pandas as pd

import bgcArgo

wmo = [ 5904852, 5904689, 5904179, 5904682, 4902414, 4901779,
        4902481, 1901879, 6902905, 6903237, 5903273, 6902805 ]
dac = 4*['aoml'] + 3*['meds'] + ['bodc'] + 4*['coriolis']

dacpath = '/ifremer/argo/dac'
fltpath = ['{}/{}/{}'.format(dacpath, d, w) for d, w in zip(dac, wmo)]

local_path = '/Users/gordonc/Documents/data/Argo'

bgcArgo.io.get_argo(fltpath, local_path=local_path)