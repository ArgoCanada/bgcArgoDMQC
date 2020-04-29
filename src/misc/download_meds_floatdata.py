#!/usr/bin/python

import sys
import ftplib

from pathlib import Path

import numpy as np
import pandas as pd

# load in MEDS list of floats with doxy sensor and metadata
p = Path('../../data/external/')
meds_doxy_floatlist = Path('../../data/external/meds_doxy_floats.csv')
# load the data into a dataframe
df = pd.read_csv(meds_doxy_floatlist)
# remove leading Q from float numbers and any spaces
df.WMO_NUMBER = np.array([n.strip('Q').strip(' ') for n in df.WMO_NUMBER])
# remove coriolis datacentre we only want meds
df = df.drop(df.index[df.PROCESSED_BY == 'Coriolis'])

# connect to dac ftp server
url = 'ftp.ifremer.fr'
dacdir = '/ifremer/argo/dac/meds'
# ftp object
ftp = ftplib.FTP(url)
ftp.login()
ftp.cwd(dacdir)

# loop through each float
for flt in df.WMO_NUMBER:
    # ------------------ SUMMARY FILES (Sprof, BRtraj, meta) ------------------
    
    # change ftp directory
    ftp.cwd(flt)
    # this is basically the same as 'ls' in unix
    files = ftp.nlst('*.nc')
    # define local location to save file
    local_path = Path('../../data/meds/') / flt
    
    # make the directory if it doesn't exist
    if not local_path.is_dir():
        local_path.mkdir()
    
    # download the files
    for fn in files:
        # define the local file to have the same name as on the FTP server
        local_file = local_path / fn
        # only download the file if it doesn't already exist locally
        if not local_file.exists():
            print(local_file)
            # open the local file
            lf = open(local_file, 'wb')
            # retrieve the file on FTP server,
            ftp.retrbinary('RETR ' + fn, lf.write)

    # ------------------------ INDIVIDUAL PROFILE FILES -----------------------
    # repeat as above
    if 'profiles' in ftp.nlst():
        ftp.cwd('profiles')
        files = ftp.nlst('*.nc')

        local_path = local_path / 'profiles'
        if not local_path.exists():
            local_path.mkdir()

        for fn in files:
            local_file = local_path / fn
            if not local_file.exists():
                print(local_file)
                lf = open(local_file, 'wb')
                ftp.retrbinary('RETR ' + fn, lf.write)

        # back to parent directory
        ftp.cwd('../../')
    else:
        ftp.cwd('../')
