#!/usr/bin/python

from pathlib import Path
import numpy as np

gfn = Path('../data/file_gains_to_be_processed.npy')
fg  = np.load(gfn, allow_pickle=True)

mfn = Path('../data/file_comments_to_be_processed.npy')
fm  = np.load(mfn, allow_pickle=True)

for glist in fg:
    if type(glist) is float:
        g = glist
        print(g)
    else:
        for gstr in glist:
            g = float(gstr.split('=')[-1].strip())
            print(g)