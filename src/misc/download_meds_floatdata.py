#!/usr/bin/python

import os
import sys
import ftplib

import numpy as np
import pandas as pd

# load in MEDS list of floats with doxy sensor and metadata
meds_doxy_floatlist = '../data/external/meds_doxy_floats.csv'
df = pd.read_csv(meds_doxy_floatlist)
# remove leading Q from float numbers and any spaces
df.WMO_NUMBER = np.array([n.strip('Q').strip(' ') for n in df.WMO_NUMBER])
df = df.drop(df.index[df.PROCESSED_BY == 'Coriolis'])

# connect to dac ftp server
url = 'ftp.ifremer.fr'
dacdir = '/ifremer/argo/dac/meds'

ftp = ftplib.FTP(url)
ftp.login()
ftp.cwd(dacdir)

for flt in df.WMO_NUMBER:
    ftp.cwd(flt)
    files = ftp.nlst('*.nc')

    local_path = os.path.join('../data/meds/',flt)
    if not os.path.exists(local_path):
        os.makedirs(local_path)

    for fn in files:
        local_file = os.path.join(local_path,fn)
        if not os.path.isfile(local_file):
            print(local_file)
            lf = open(local_file, 'wb')
            ftp.retrbinary('RETR ' + fn, lf.write)

    if 'profiles' in ftp.nlst():
        ftp.cwd('profiles')
        files = ftp.nlst('*.nc')

        local_path = os.path.join(local_path, 'profiles')
        if not os.path.exists(local_path):
            os.makedirs(local_path)

        for fn in files:
            local_file = os.path.join(local_path,fn)
            if not os.path.isfile(local_file):
                print(local_file)
                lf = open(local_file, 'wb')
                ftp.retrbinary('RETR ' + fn, lf.write)

        ftp.cwd('../../')
    else:
        ftp.cwd('../')
