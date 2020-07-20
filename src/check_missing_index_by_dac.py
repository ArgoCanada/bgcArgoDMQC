#!/usr/bin/python

import sys
import bgcArgo as bgc

index = bgc.get_index(index='global')
# for dates or location
for dac in index.dac.unique():
    sub_ix = index[index.dac == dac]
    missing_ix = sub_ix.date.isnull() | sub_ix.latitude.isnull() | sub_ix.longitude.isnull()
    N_missing = missing_ix.sum()
    sys.stdout.write('{:.5f}%\t{}\n'.format(100*N_missing/sub_ix.shape[0], dac))

# for just dates
for dac in index.dac.unique():
    sub_ix = index[index.dac == dac]
    missing_ix = sub_ix.date.isnull()
    N_missing = missing_ix.sum()
    sys.stdout.write('{:.5f}%\t{}\n'.format(100*N_missing/sub_ix.shape[0], dac))

# get the meds one
sub_ix = index[index.dac == 'meds']
missing_ix = sub_ix.date.isnull()
df = sub_ix[missing_ix]
