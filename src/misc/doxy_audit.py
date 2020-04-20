#!/usr/bin/python

import os

import numpy as np
import pandas as pd

audit_file = '../data/external/DOXY_audit_022020.TXT'
gain_file  = '../data/external/DOXY_gains_022020.TXT'

af = pd.read_csv(audit_file, sep='\t', header=20)
# gf = pd.read_csv(gain_file, sep='\t', header=19)

# full audit from February
# get only floats with MEDS as their DAC
meds = af.drop(af.index[af.DAC != 'meds'])
# how many different floats are in audit?
audit_floats = meds.WMO.unique()

# meds-only audit file from April
new_audit_file = '../data/external/DOXY_audit_meds_041620.TXT'
df = pd.read_csv(new_audit_file, sep='\t', header=20)
