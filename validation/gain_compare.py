#!/usr/bin/python

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from bgcArgo import sprof, profiles

sns.set(context='paper', style='whitegrid', palette='colorblind')
sage_path = Path('/Users/gordonc/Documents/data/Argo/sage/')
sage_files = list(sage_path.glob('*float*.h5'))

fn = sage_files[0]
wmo = int(str(fn).split('\\')[-1].split('_')[-1].split('.')[0])

datapath = Path('/Users/gordonc/Documents/data/')
profiles.set_dirs(argo_path=datapath / 'Argo', woa_path=datapath / 'WOA18', ncep_path=datapath / 'NCEP')
syn = sprof(wmo)
gains = syn.calc_gains(ref='WOA')