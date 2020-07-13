#!/usr/bin/python

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import bgcArgo as bgc

# load in DOXY audit file - most recent one was on July 9 2020
audit_file = Path('../data/DOXY_audit_070720.TXT')
df = pd.read_csv(audit_file, sep='\t', header=25)

# download the synthetic, meta, and BRtraj files for each float in the audit
dacpath = '/ifremer/argo/dac'
fltpath = ['{}/{}/{}'.format(dacpath, bgc.get_dac(w), w) for w in df.WMO.unique()]
local_path = '/Users/gordonc/Documents/data/Argo'
bgc.io.get_argo(fltpath, local_path=local_path, mode='summary')
