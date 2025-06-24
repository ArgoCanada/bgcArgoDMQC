import sys
import warnings

from pathlib import Path
import ftplib

from time import time
import numpy as np
import pandas as pd

from .. import resource

def index_exists():

    meta  = 'ar_index_global_meta.txt.gz'
    index = 'ar_index_global_prof.txt.gz'
    bgc   = 'argo_bio-profile_index.txt.gz'
    synth = 'argo_synthetic-profile_index.txt.gz'
    traj  = 'ar_index_global_traj.txt.gz'
    local_meta  = resource.path('Index') / meta
    local_index = resource.path('Index') / index
    local_bgc   = resource.path('Index') / bgc
    local_synth = resource.path('Index') / synth
    local_traj  = resource.path('Index') / traj

    return all([local_meta.exists(), local_index.exists(), local_bgc.exists(), local_synth.exists(), local_traj.exists()])

def read_index(mission='B', source='argopy', url=resource.URL):
    '''
    Function to read and extract information from Argo global index,
    then save it to a dataframe for faster access.

    Args:
        mission (str): *bgc-b*, *core*, *bgc-s*, *traj*, or *meta* for biogeochemical, global/core, synthetic, and metadata indices respectfully. 
    '''

    if source != 'argopy':
        url_dir = resource.URL_DIR_DICT[url]
        if mission == 'bgc-b':
            local_filename = resource.path('Index') / 'argo_bio-profile_index.txt.gz'
            remote_filename = 'ftp://' + url + (Path(url_dir) / 'argo_bio-profile_index.txt.gz').as_posix()
        elif mission == 'core':
            local_filename = resource.path('Index') / 'ar_index_global_prof.txt.gz'
            remote_filename = 'ftp://' + url + (Path(url_dir) / 'ar_index_global_prof.txt.gz').as_posix()
        elif mission == 'bgc-s':
            local_filename = resource.path('Index') / 'argo_synthetic-profile_index.txt.gz'
            remote_filename = 'ftp://' + url + (Path(url_dir) / 'argo_synthetic-profile_index.txt.gz').as_posix()
        elif mission == 'meta':
            local_filename = resource.path('Index') / 'ar_index_global_meta.txt.gz'
            remote_filename = 'ftp://' + url + (Path(url_dir) / 'ar_index_global_meta.txt.gz').as_posix()
        elif mission == 'traj':
            local_filename = resource.path('Index') / 'ar_index_global_traj.txt.gz'
            remote_filename = 'ftp://' + url + (Path(url_dir) / 'ar_index_global_traj.txt.gz').as_posix()
        else: # pragma: no cover
            raise ValueError(f'Input {mission} not recognized - must be *bgc-b*, *core*, *bgc-s*, *traj*, or *meta*')

    if source == 'argopy':
        if mission == 'traj' or mission == 'meta':
            raise ValueError('argopy does not support traj or meta file handling yet - use source="local" or source="remote" to load.')
        import argopy
        df = argopy.ArgoIndex(index_file=mission).to_dataframe()
    elif source == 'remote':
        df = pd.read_csv(remote_filename, compression='gzip', header=8)

        df['wmo'] = np.array([int(f.split('/')[1]) for f in df.file])
        df['dac'] = np.array([f.split('/')[0] for f in df.file])
        if mission != 'meta' and mission != 'traj':
            df['cycle'] = np.array([int(f.split('/')[-1].split('.')[-2].split('_')[-1].replace('D','')) for f in df.file])

    elif source == 'local':
        if not Path(local_filename).exists():
            sys.stdout.write('Index file does not exist, downloading now, this may take a few minutes\n')
            update_index(ftype=mission)

        df =  pd.read_csv(local_filename, compression='gzip', header=8)

        df['wmo'] = np.array([int(f.split('/')[1]) for f in df.file])
        df['dac'] = np.array([f.split('/')[0] for f in df.file])
        if mission != 'meta' and mission != 'traj':
            df['cycle'] = np.array([int(f.split('/')[-1].split('.')[-2].split('_')[-1].replace('D','')) for f in df.file])

    return df

def update_index(ftype=None, url=resource.URL):
    '''
    Function to access FTP server to download Argo metadata and profile global
    index files
    '''

    ftp = ftplib.FTP(url)
    url_dir = resource.URL_DIR_DICT[url]
    ftp.login()
    ftp.cwd(url_dir)

    meta  = 'ar_index_global_meta.txt.gz'
    traj  = 'ar_index_global_traj.txt.gz'
    index = 'ar_index_global_prof.txt.gz'
    bgc   = 'argo_bio-profile_index.txt.gz'
    synth = 'argo_synthetic-profile_index.txt.gz'


    local_meta = resource.path('Index') / meta
    if ftype is None or ftype == 'meta':
        lf = open(local_meta, 'wb')
        ftp.retrbinary('RETR ' + meta, lf.write)
        lf.close()

    local_traj = resource.path('Index') / traj
    if ftype is None or ftype == 'traj':
        lf = open(local_traj, 'wb')
        ftp.retrbinary('RETR ' + traj, lf.write)
        lf.close()

    local_index = resource.path('Index') / index
    if ftype is None or ftype =='profile' or ftype == 'C':
        lf = open(local_index, 'wb')
        ftp.retrbinary('RETR ' + index, lf.write)
        lf.close()

    local_bgc = resource.path('Index') / bgc
    if ftype is None or ftype =='bgc' or ftype == 'B':
        lf = open(local_bgc, 'wb')
        ftp.retrbinary('RETR ' + bgc, lf.write)
        lf.close()

    local_synth = resource.path('Index') / synth
    if ftype is None or ftype =='synthetic' or ftype == 'S':
        lf = open(local_synth, 'wb')
        ftp.retrbinary('RETR ' + synth, lf.write)
        lf.close()

    return ftp

def check_index_time(local_meta, local_index, local_bgc, local_synth, local_traj):
    meta_mtime  = local_meta.stat().st_mtime
    index_mtime = local_index.stat().st_mtime
    bgc_mtime   = local_bgc.stat().st_mtime
    synth_mtime = local_synth.stat().st_mtime
    traj_mtime  = local_traj.stat().st_mtime

    # get time since last update
    curr_time   = time()
    curr_time   = time()
    meta_delta  = curr_time - meta_mtime
    index_delta = curr_time - index_mtime
    bgc_delta   = curr_time - bgc_mtime
    synth_delta = curr_time - synth_mtime
    traj_delta = curr_time - traj_mtime

    if meta_delta / 60 / 60 / 24 > 7:
        d = meta_delta / 60 / 60 / 24
        sys.stdout.write('Argo global metadata index is more than 7 days old - has not been updated in {:d} days  - downloading now - this may take some time depending on your internet connection\n'.format(int(d)))
        update_index(ftype='meta')

    if index_delta / 60 / 60 / 24 > 7:
        d = index_delta / 60 / 60 / 24
        sys.stdout.write('Argo global profile index is more than 7 days old - has not been updated in {:d} days  - downloading now - this may take some time depending on your internet connection\n'.format(int(d)))
        update_index(ftype='profile')

    if bgc_delta / 60 / 60 / 24 > 7:
        d = bgc_delta / 60 / 60 / 24
        sys.stdout.write('Argo global BGC index is more than 7 days old - has not been updated in {:d} days  - downloading now - this may take some time depending on your internet connection\n'.format(int(d)))
        update_index(ftype='bgc')

    if synth_delta / 60 / 60 / 24 > 7:
        d = synth_delta / 60 / 60 / 24
        sys.stdout.write('Argo global synthetic profile index is more than 7 days old - has not been updates in {:d} days  - downloading now - this may take some time depending on your internet connection\n'.format(int(d)))
        update_index(ftype='synthetic')

    if traj_delta / 60 / 60 / 24 > 7:
        d = synth_delta / 60 / 60 / 24
        sys.stdout.write('Argo global synthetic profile index is more than 7 days old - has not been updates in {:d} days  - downloading now - this may take some time depending on your internet connection\n'.format(int(d)))
        update_index(ftype='synthetic')

def check_index(mode=None):
    '''
    Function to check age of Argo metadata and profile global index files,
    warn if they have not been updated in more than 1 week. Runs on import.
    '''

    if mode is None:
        meta  = 'ar_index_global_meta.txt.gz'
        index = 'ar_index_global_prof.txt.gz'
        bgc   = 'argo_bio-profile_index.txt.gz'
        synth = 'argo_synthetic-profile_index.txt.gz'
        traj  = 'ar_index_global_traj.txt.gz'
        local_meta  = resource.path('Index') / meta
        local_index = resource.path('Index') / index
        local_bgc   = resource.path('Index') / bgc
        local_synth = resource.path('Index') / synth
        local_traj  = resource.path('Index') / traj

        check_index_time(local_meta, local_index, local_bgc, local_synth, local_traj)

    elif mode == 'install':
        meta  = 'ar_index_global_meta.txt.gz'
        index = 'ar_index_global_prof.txt.gz'
        bgc   = 'argo_bio-profile_index.txt.gz'
        synth = 'argo_synthetic-profile_index.txt.gz'
        traj  = 'ar_index_global_traj.txt.gz'
        local_meta  = resource.path('Index') / meta
        local_index = resource.path('Index') / index
        local_bgc   = resource.path('Index') / bgc
        local_synth = resource.path('Index') / synth
        local_traj  = resource.path('Index') / traj

        if not ((local_meta.exists() and local_index.exists()) and (local_bgc.exists() and local_synth.exists())):
            sys.stdout.write('At least one index file does not exist - downloading now - this may take some time depending on your internet connection\n')
            update_index()
        else: # pragma: no cover
            check_index_time(local_meta, local_index, local_bgc, local_synth, local_traj)

def get_dac(wmo):

    if '__globalindex__' not in globals():
            global __globalindex__
            __globalindex__ = read_index(mission='C')
    
    dac = __globalindex__.loc[__globalindex__.wmo == wmo].dac.iloc[0]

    return dac