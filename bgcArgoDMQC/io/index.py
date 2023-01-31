import sys
import warnings

from pathlib import Path
import ftplib

from time import time
import numpy as np
import pandas as pd

from .. import configure

global index_path
index_path = Path(__file__).parent.parent.absolute() / 'ref'
index_path.mkdir(exist_ok=True)

global URL_DICT
URL_DICT = {
    'ftp.ifremer.fr':'ftp.ifremer.fr', 
    'ifremer':'ftp.ifremer.fr', 
    'coriolis':'ftp.ifremer.fr', 
    'usgodae.org':'usgodae.org', 
    'godae':'usgodae.org', 
    'us':'usgodae.org'
}

global URL_DIR_DICT
URL_DIR_DICT = {
    'ftp.ifremer.fr':'/ifremer/argo/', 
    'usgodae.org':'/pub/outgoing/argo/', 
}

global URL
config = configure.read_config()
if 'default_url' in config.keys():
    url_name = config.pop('default_url')
    URL = URL_DICT[url_name]
else:
    URL = 'ftp.ifremer.fr'

def index_exists():

    meta  = 'ar_index_global_meta.txt.gz'
    index = 'ar_index_global_prof.txt.gz'
    bgc   = 'argo_bio-profile_index.txt.gz'
    synth = 'argo_synthetic-profile_index.txt.gz'
    local_meta  = index_path / meta
    local_index = index_path / index
    local_bgc   = index_path / bgc
    local_synth = index_path / synth

    return all([local_meta.exists(), local_index.exists(), local_bgc.exists(), local_synth.exists()])

def read_index(mission='B', remote=False, url=URL):
    '''
    Function to read and extract information from Argo global index,
    then save it to a dataframe for faster access.

    Args:
        mission (str): *B*, *C*, *S*, or *M* for biogeochemical, global/core, synthetic, and metadata indices respectfully. 
    '''

    url_dir = URL_DIR_DICT[url]
    if mission == 'B':
        local_filename = index_path / 'argo_bio-profile_index.txt.gz'
        remote_filename = 'ftp://' + url + (Path(url_dir) / 'argo_bio-profile_index.txt.gz').as_posix()
    elif mission == 'C':
        local_filename = index_path / 'ar_index_global_prof.txt.gz'
        remote_filename = 'ftp://' + url + (Path(url_dir) / 'ar_index_global_prof.txt.gz').as_posix()
    elif mission == 'S':
        local_filename = index_path / 'argo_synthetic-profile_index.txt.gz'
        remote_filename = 'ftp://' + url + (Path(url_dir) / 'argo_synthetic-profile_index.txt.gz').as_posix()
    elif mission == 'M':
        local_filename = index_path / 'ar_index_global_meta.txt.gz'
        remote_filename = 'ftp://' + url + (Path(url_dir) / 'ar_index_global_meta.txt.gz').as_posix()
    elif mission == 'T':
        local_filename = index_path / 'ar_index_global_traj.txt.gz'
        remote_filename = 'ftp://' + url + (Path(url_dir) / 'ar_index_global_traj.txt.gz').as_posix()
    else:
        raise ValueError('Input {} not recognized'.format(mission))

    if remote:
        df = pd.read_csv(remote_filename, compression='gzip', header=8)
    else:
        # warnings.warn('Could not read index file from ifremer FTP server, trying to load using local file, which may not be up to date.')
        if not Path(local_filename).exists():
            sys.stdout.write('Index file does not exist, downloading now, this may take a few minutes\n')
            update_index(ftype=mission)

        df =  pd.read_csv(local_filename, compression='gzip', header=8)

    df['wmo'] = np.array([int(f.split('/')[1]) for f in df.file])
    df['dac'] = np.array([f.split('/')[0] for f in df.file])
    if mission != 'M' and mission != 'T':
        df['cycle'] = np.array([int(f.split('/')[-1].split('.')[-2].split('_')[-1].replace('D','')) for f in df.file])

    return df

def update_index(ftype=None, url=URL):
    '''
    Function to access FTP server to download Argo metadata and profile global
    index files
    '''

    ftp = ftplib.FTP(url)
    url_dir = URL_DIR_DICT[url]
    ftp.login()
    ftp.cwd(url_dir)

    meta  = 'ar_index_global_meta.txt.gz'
    traj  = 'ar_index_global_traj.txt.gz'
    index = 'ar_index_global_prof.txt.gz'
    bgc   = 'argo_bio-profile_index.txt.gz'
    synth = 'argo_synthetic-profile_index.txt.gz'


    local_meta = index_path / meta
    if ftype is None or ftype == 'meta':
        lf = open(local_meta, 'wb')
        ftp.retrbinary('RETR ' + meta, lf.write)
        lf.close()

    local_traj = index_path / traj
    if ftype is None or ftype == 'traj':
        lf = open(local_traj, 'wb')
        ftp.retrbinary('RETR ' + traj, lf.write)
        lf.close()

    local_index = index_path / index
    if ftype is None or ftype =='profile' or ftype == 'C':
        lf = open(local_index, 'wb')
        ftp.retrbinary('RETR ' + index, lf.write)
        lf.close()

    local_bgc = index_path / bgc
    if ftype is None or ftype =='bgc' or ftype == 'B':
        lf = open(local_bgc, 'wb')
        ftp.retrbinary('RETR ' + bgc, lf.write)
        lf.close()

    local_synth = index_path / synth
    if ftype is None or ftype =='synthetic' or ftype == 'S':
        lf = open(local_synth, 'wb')
        ftp.retrbinary('RETR ' + synth, lf.write)
        lf.close()

    return ftp

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
        local_meta  = index_path / meta
        local_index = index_path / index
        local_bgc   = index_path / bgc
        local_synth = index_path / synth

        meta_mtime  = local_meta.stat().st_mtime
        index_mtime = local_index.stat().st_mtime
        bgc_mtime   = local_bgc.stat().st_mtime
        synth_mtime = local_index.stat().st_mtime

        # get time since last update
        curr_time   = time()
        meta_delta  = curr_time - meta_mtime
        index_delta = curr_time - index_mtime
        bgc_delta   = curr_time - bgc_mtime
        synth_delta = curr_time - synth_mtime

        if meta_delta / 60 / 60 / 24 > 7:
            d = meta_delta / 60 / 60 / 24
            warnings.warn('Argo global metadata index is more than 7 days old - has not been updated in {:d} days - consider running update_index()'.format(int(d)), Warning)

        if index_delta / 60 / 60 / 24 > 7:
            d = index_delta / 60 / 60 / 24
            warnings.warn('Argo global profile index is more than 7 days old - has not been updated in {:d} days - consider running update_index()'.format(int(d)), Warning)

        if bgc_delta / 60 / 60 / 24 > 7:
            d = bgc_delta / 60 / 60 / 24
            warnings.warn('Argo global BGC index is more than 7 days old - has not been updated in {:d} days - consider running update_index()'.format(int(d)), Warning)

        if synth_delta / 60 / 60 / 24 > 7:
            d = synth_delta / 60 / 60 / 24
            warnings.warn('Argo global synthetic profile index is more than 7 days old - has not been updated in {:d} days - consider running update_index()'.format(int(d)), Warning)

    elif mode == 'install':
        meta  = 'ar_index_global_meta.txt.gz'
        index = 'ar_index_global_prof.txt.gz'
        bgc   = 'argo_bio-profile_index.txt.gz'
        synth = 'argo_synthetic-profile_index.txt.gz'
        local_meta  = index_path / meta
        local_index = index_path / index
        local_bgc   = index_path / bgc
        local_synth = index_path / synth

        if not ((local_meta.exists() and local_index.exists()) and (local_bgc.exists() and local_synth.exists())):
            sys.stdout.write('At least one index file does not exist - downloading now - this may take some time depending on your internet connection\n')
            update_index()
        else:
            meta_mtime  = local_meta.stat().st_mtime
            index_mtime = local_index.stat().st_mtime
            bgc_mtime   = local_bgc.stat().st_mtime
            synth_mtime = local_synth.stat().st_mtime

            # get time since last update
            curr_time   = time()
            curr_time   = time()
            meta_delta  = curr_time - meta_mtime
            index_delta = curr_time - index_mtime
            bgc_delta   = curr_time - bgc_mtime
            synth_delta = curr_time - synth_mtime

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

def get_dac(wmo):

    if '__globalindex__' not in globals():
            global __globalindex__
            __globalindex__ = read_index(mission='C')
    
    dac = __globalindex__.loc[__globalindex__.wmo == wmo].dac.iloc[0]

    return dac