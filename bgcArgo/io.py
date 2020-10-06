#!/usr/bin/python

import sys
import warnings

from time import time
from pathlib import Path
import ftplib
import pandas as pd

import numpy as np
import matplotlib.dates as mdates

from netCDF4 import Dataset

from . import util

global index_path
index_path = Path(__file__).parent.absolute() / 'ref'
if not index_path.exists():
    index_path.mkdir()

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

def read_index(mission='B'):
    '''
    Function to read and extract information from Argo global index,
    then save it to a dataframe for faster access
    '''

    if mission == 'B':
        filename = index_path / 'argo_bio-profile_index.txt.gz'
    elif mission == 'C':
        filename = index_path / 'ar_index_global_prof.txt.gz'
    elif mission == 'S':
        filename = index_path / 'argo_synthetic-profile_index.txt.gz'

    if not Path(filename).exists():
        sys.stdout.write('Index file does not exist, downloading now, this may take a few minutes\n')
        update_index(ftype=mission)

    df =  pd.read_csv(filename, compression='gzip', header=8)

    df['dac'] = np.array([f.split('/')[0] for f in df.file])
    df['wmo'] = np.array([int(f.split('/')[1]) for f in df.file])
    df['cycle'] = np.array([int(f.split('/')[-1].split('.')[-2].split('_')[-1].replace('D','')) for f in df.file])

    return df

def update_index(ftype=None):
    '''
    Function to access FTP server to download Argo metadata and profile global
    index files
    '''

    ftp = ftplib.FTP('ftp.ifremer.fr')
    ftp.login()
    ftp.cwd('/ifremer/argo/')

    meta  = 'ar_index_global_meta.txt.gz'
    index = 'ar_index_global_prof.txt.gz'
    bgc   = 'argo_bio-profile_index.txt.gz'
    synth = 'argo_synthetic-profile_index.txt.gz'

    local_meta = index_path / meta
    if ftype is None or ftype == 'meta':
        lf = open(local_meta, 'wb')
        ftp.retrbinary('RETR ' + meta, lf.write)
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
    
    dac = __globalindex__[__globalindex__.wmo == wmo].dac.iloc[0]

    return dac

def get_woa18(varname, local_path='./', ftype='netcdf', overwrite=False):
    '''
    Function to download WOA data for a given variable

    INPUT:
              varname: woa18 variable to download data for - one of:
                  T: temperature
                  S: salinity
                  O2: dissolved oxygen
                  O2sat: oxygen percent saturation
                  NO3: nitrate
                  Si: silicate
                  PO4: phosphate
              local_path: path to save files to (inside variable folder, ex.
                  if local_path = '../data/woa', and varname = 'T', files
                  will be saved in ../data/woa/temperature/*.nc. Defaults
                  to current directory
              ftype: file format to download, defaults to netcdf (.nc):
                  ascii: .dat
                  csv: .csv
                  netcdf: .nc
              overwrite: boolean flag, if False, does not re-download
                  existing files, if true, will download no matter what,
                  defaults to False

    OUTPUT:
              ftp: the ftplib server object

    AUTHOR:   Christopher Gordon
              Fisheries and Oceans Canada
              chris.gordon@dfo-mpo.gc.ca

    LAST UPDATE: 29-04-2020

    CHANGE LOG:
    '''

    local_path = Path(local_path)
    url = 'ftp.nodc.noaa.gov'
    param, dtype, ftpdir = util.decode_woa_var(varname)

    ftp = ftplib.FTP(url, 'anonymous', 'chris.gordon@dfo-mpo.gc.ca')
    ftp.cwd('pub/woa/WOA18/DATA/{}/{}/{}/1.00/'.format(ftpdir, ftype, dtype))

    local_path = local_path / ftpdir

    if not local_path.is_dir():
        local_path.mkdir()

    for fn in ftp.nlst():
        local_file = local_path / fn
        if not local_file.exists() or overwrite:
            print(local_file)
            # open the local file
            lf = open(local_file, 'wb')
            # retrieve the file on FTP server,
            ftp.retrbinary('RETR ' + fn, lf.write)
            lf.close()

    return ftp

def get_ncep(varname, local_path='./', overwrite=False, years=[2010, 2020]):
    '''
    Function to download NCEP reanalysis gaussian gridded surface air
    pressure data

    INPUT:
              varname: 'pres' (pressure) or 'rhum' (relative humidity)
                  or 'land' (to get land mask)
              local_path: path to save files to, defaults
                  to current directory
              overwrite: boolean flag, if False, does not re-download
                  existing files, if true, will download no matter what,
                  defaults to False

    OUTPUT:
              ftp: the ftplib server object

    AUTHOR:   Christopher Gordon
              Fisheries and Oceans Canada
              chris.gordon@dfo-mpo.gc.ca

    LAST UPDATE: 29-04-2020

    CHANGE LOG:
    '''

    local_path = Path(local_path)
    url = 'ftp.cdc.noaa.gov'

    ftp = ftplib.FTP(url, 'anonymous')

    if type(years) is int:
        yearlist = [years]
    elif len(years) == 1:
        yearlist = years
    elif len(years) == 2:
        yearlist = range(years[0], years[1]+1)
    elif len(years > 2):
        yearlist = years
    
    if varname == 'pres':
        ftp.cwd('Datasets/ncep.reanalysis2/gaussian_grid/')

        local_path = local_path / varname
        if not local_path.is_dir():
            local_path.mkdir()

        for yr in yearlist:
            fn = 'pres.sfc.gauss.{}.nc'.format(yr)
            local_file = local_path / fn

            if not local_file.exists() or overwrite:
                print(local_file)
                # open the local file
                lf = open(local_file, 'wb')
                # retrieve the file on FTP server,
                ftp.retrbinary('RETR ' + fn, lf.write)
                lf.close()

    elif varname == 'rhum':

        ftp.cwd('Datasets/ncep.reanalysis/surface/')

        local_path = local_path / varname
        if not local_path.is_dir():
            local_path.mkdir()

        for yr in yearlist:
            fn = 'rhum.sig995.{}.nc'.format(yr)
            local_file = local_path / fn
            if not local_file.exists() or overwrite:
                print(local_file)
                # open the local file
                lf = open(local_file, 'wb')
                # retrieve the file on FTP server,
                ftp.retrbinary('RETR ' + fn, lf.write)
                lf.close()

    elif varname == 'land':

        local_path = local_path / varname
        if not local_path.is_dir():
            local_path.mkdir()

        ftp.cwd('Datasets/ncep.reanalysis2/gaussian_grid/')
        fn = 'land.sfc.gauss.nc'
        local_file = local_path / fn
        if not local_file.exists() or overwrite:
            lf = open(local_file, 'wb')
            ftp.retrbinary('RETR ' + fn, lf.write)
            lf.close()

        ftp.cwd('../../ncep.reanalysis/surface/')
        fn = 'land.nc'
        local_file = local_path / fn
        if not local_file.exists() or overwrite:
            lf = open(local_file, 'wb')
            ftp.retrbinary('RETR ' + fn, lf.write)
            lf.close()

    else:
        raise ValueError('Invalid varname input')

    return ftp

def get_argo(*args, local_path='./', url='ftp.ifremer.fr', overwrite=False, ftype=None, mission=None, mode='RD'):
    '''
    Function to download all data from a single float, or individual
    profiles

    INPUT:
              Inputs may vary depending on desired performance. Multiple
              arguments may be provided to download all files from a certain
              float or argo defined geographical area. A single path to a
              file may be provided to download that file. A list of files
              may be provided as well.

              overwrite: boolean flag, if False, does not re-download
                  existing files, if true, will download no matter what,
                  defaults to False

    OUTPUT:

    AUTHOR:   Christopher Gordon
              Fisheries and Oceans Canada
              chris.gordon@dfo-mpo.gc.ca

    LAST UPDATE: 29-04-2020

    CHANGE LOG:
    '''

    if len(args) == 1:
        local_path = Path(local_path)

        ftp = ftplib.FTP(url)
        ftp.login()

        if type(args[0]) is not list:
            args[0] = list(args[0])

        for init_a in args[0]:
            # if its a float number, build ftp paths to floats
            if type(init_a) is int or type(init_a) is float:
                if url == 'ftp.ifremer.fr':
                    base_path = '/ifremer/argo/dac'
                elif url == 'usgodae.org':
                    base_path = '/pub/outgoing/argo/dac'

                a = '{}/{}/{:d}'.format(base_path, get_dac(init_a), init_a)
            else:
                a = init_a

            # if its a file
            if a[-2:] == 'nc':
                # if its a profile file
                if 'profiles' in a:
                    ftp_wmo_path = '/'.join(a.split('/')[:-2])
                    wmo = ftp_wmo_path.split('/')[-1]
                    dac = ftp_wmo_path.split('/')[-2]
                    ftp.cwd(ftp_wmo_path)
                    ftp.cwd('profiles')
                    fn = a.split('/')[-1]

                    # define local location to save file
                    dac_path = local_path / dac
                    wmo_path = local_path / dac / wmo
                    profile_path = wmo_path / 'profiles'

                    # make the directory if it doesn't exist
                    if not dac_path.is_dir():
                        dac_path.mkdir()
                    if not wmo_path.is_dir():
                        wmo_path.mkdir()
                    if not profile_path.is_dir():
                        profile_path.mkdir()

                    # define the local file to have the same name as on the FTP server
                    wmo_file = profile_path / fn
                    # only download the file if it doesn't already exist locally
                    if not wmo_file.exists() or overwrite:
                        print(wmo_file)
                        # open the local file
                        lf = open(wmo_file, 'wb')
                        # retrieve the file on FTP server,
                        ftp.retrbinary('RETR ' + fn, lf.write)
                        lf.close()

                # if not - so Sprof, Mprof, BRtraj or meta
                else:
                    ftp_wmo_path = '/'.join(a.split('/')[:-1])
                    wmo = ftp_wmo_path.split('/')[-1]
                    dac = ftp_wmo_path.split('/')[-2]
                    ftp.cwd(ftp_wmo_path)
            # if its a float directory
            else:
                ftp_wmo_path = a
                wmo = ftp_wmo_path.split('/')[-1]
                dac = ftp_wmo_path.split('/')[-2]
                ftp.cwd(ftp_wmo_path)

                files = ftp.nlst('*.nc')
                    
                # define local location to save file
                dac_path = local_path / dac
                wmo_path = local_path / dac / wmo

                # make the directory if it doesn't exist
                if not dac_path.is_dir():
                    dac_path.mkdir()
                if not wmo_path.is_dir():
                    wmo_path.mkdir()
                

                # download the files
                for fn in files:
                    # define the local file to have the same name as on the FTP server
                    wmo_file = wmo_path / fn
                    # only download the file if it doesn't already exist locally
                    if not wmo_file.exists() or overwrite:
                        print(wmo_file)
                        # open the local file
                        lf = open(wmo_file, 'wb')
                        # retrieve the file on FTP server,
                        ftp.retrbinary('RETR ' + fn, lf.write)
                        lf.close()

                # repeat as above
                if 'profiles' in ftp.nlst() and ftype != 'summary':
                    ftp.cwd('profiles')
                    if mission is None or mission == 'CB':
                        if mode == 'RD':
                            files = ftp.nlst('*.nc')
                        elif mode == 'R':
                            files = ftp.nlst('*.nc')
                            ix = np.array(['D' in fn for fn in files])
                            files = list(np.array(files)[~ix])
                        elif mode == 'D':
                            files = ftp.nlst('D*.nc') + ftp.nlst('BD*.nc')
                    if mission == 'C':
                        if mode == 'RD':
                            files = ftp.nlst('*.nc')
                            ix = np.array(['B' in fn for fn in files])
                            files = list(np.array(files)[~ix])
                        elif mode == 'R':
                            files = ftp.nlst('*.nc')
                            ix = np.array(['B' in fn or 'D' in fn for fn in files])
                            files = list(np.array(files)[~ix])
                        elif mode == 'D':
                            files = ftp.nlst('D*.nc')
                    if mission == 'B':
                        if mode == 'RD':
                            files = ftp.nlst('B*.nc')
                        elif mode == 'R':
                            files = ftp.nlst('B*.nc')
                            ix = np.array(['D' in fn for fn in files])
                            files = list(np.array(files)[~ix])
                        elif mode == 'D':
                            files = ftp.nlst('BD*.nc')

                    profile_path = wmo_path / 'profiles'
                    if not profile_path.exists():
                        profile_path.mkdir()

                    for fn in files:
                        profile_file = profile_path / fn
                        if not profile_file.exists() or overwrite:
                            print(profile_file)
                            lf = open(profile_file, 'wb')
                            ftp.retrbinary('RETR ' + fn, lf.write)
                            lf.close()

    elif len(args) > 1:
        local_path = Path(local_path)

        dacdir = args[0]
        ftp = ftplib.FTP(url)
        ftp.login()
        ftp.cwd(dacdir)

        wmo_numbers = args[1]
        for flt in wmo_numbers:
            if type(flt) is str:
                wmo = flt
            else:
                wmo = str(flt)

            # ------------------ SUMMARY FILES (Sprof, BRtraj, meta) --------------
            # change ftp directory
            print(wmo)
            ftp.cwd(wmo)
            files = ftp.nlst('*.nc')
            # define local location to save file
            wmo_path = local_path / wmo

            # make the directory if it doesn't exist
            if not wmo_path.is_dir():
                wmo_path.mkdir()

            # download the files
            for fn in files:
                # define the local file to have the same name as on the FTP server
                wmo_file = wmo_path / fn
                # only download the file if it doesn't already exist locally
                if not wmo_file.exists() or overwrite:
                    print(wmo_file)
                    # open the local file
                    lf = open(wmo_file, 'wb')
                    # retrieve the file on FTP server,
                    ftp.retrbinary('RETR ' + fn, lf.write)
                    lf.close()

            # ------------------------ INDIVIDUAL PROFILE FILES -------------------
            # repeat as above
            if 'profiles' in ftp.nlst() and ftype != 'summary':
                ftp.cwd('profiles')
                files = ftp.nlst('*.nc')

                profile_path = wmo_path / 'profiles'
                if not profile_path.exists():
                    profile_path.mkdir()

                for fn in files:
                    profile_file = profile_path / fn
                    if not profile_file.exists() or overwrite:
                        print(profile_file)
                        lf = open(profile_file, 'wb')
                        ftp.retrbinary('RETR ' + fn, lf.write)
                        lf.close()

                # back to parent directory
                ftp.cwd('../../')
            else:
                ftp.cwd('../')


    return ftp

def load_woa_data(track, param, zlim=(0,1000), local_path='./', verbose=True):
    '''
    Function to load WOA18 climatological data for comparison with autonomous
    floats. Data to be interpolated along the provided track (t, lat, lon).

    INPUT:
              track: array with the columns (SDN, lat, lon)
              param: requested variable, valid inputs are
                  T: temperature
                  S: salinity
                  O2: dissolved oxygen
                  O2sat: oxygen percent saturation
                  NO3: nitrate
                  Si: silicate
                  PO4: phosphate
              zlim: depth bounds (upper, lower), default to (0, 1000)
              local_path: local directory where WOA files are stored, assumes
                          current directory if no input

    OUTPUT:
              xtrack: same as track input, but adjusted lon if the track
                      crosses the 180/-180 meridian
              woa_track: list with z, lat, and lon arrays of WOA data
              data: gridded array of the input variable (month, z, lat, lon)

    AUTHOR:   Christopher Gordon
              Fisheries and Oceans Canada
              chris.gordon@dfo-mpo.gc.ca

    ACKNOWLEDGEMENT: this code is adapted from the SOCCOM SAGE_O2Argo matlab
    code, available via https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING,
    written by Tanya Maurer & Josh Plant

    LAST UPDATE: 29-04-2020

    CHANGE LOG:

    23-04-2020: changed zlim to optional input argument
    29-04-2020: switched file/path handling from os module to pathlib
    '''

    # make local_path a Path() object from a string, account for windows path
    local_path = Path(local_path)

    # check if float track crosses -180/180 meridian
    cross180 = False
    if np.max(np.abs(np.diff(track[:,2]))) > 340:
        cross180 = True
        lix = track[:,2] < 0
        lon_bounds = (np.nanmax(track[lix,2]), np.nanmin(track[~lix,2]))
    else:
        lon_bounds = (np.nanmin(track[:,2]), np.nanmax(track[:,2]))
    lat_bounds = (np.nanmin(track[:,1]), np.nanmax(track[:,1]))

    # set up extraction files, variables
    woa_param, woa_ftype, woa_dir = util.decode_woa_var(param)
    var_name  = woa_param + '_an'

    base_woa_file = 'woa18_{}_{}'.format(woa_ftype, woa_param)
    # loop through months
    for i in range(12):
        mo = i+1
        woa_file = base_woa_file + '{:02d}_01.nc'.format(mo)
        nc = Dataset(local_path / woa_dir / woa_file, 'r')

        if i == 0:
            z   = nc.variables['depth'][:]
            lat = nc.variables['lat'][:]
            lon = nc.variables['lon'][:]

            # depth boundaries
            z_ix = np.logical_and(z >= zlim[0], z <= zlim[1])
            z_ix = np.where(z_ix)[0]

            if zlim[1] > z[-1]:
                warnings.warn('Max requested depth {} greater than WOA max depth {}\n'.format(zlim[1], z[-1]), Warning)

            lat_ix = util.get_lat_index(lat, lat_bounds)
            lon_ix = util.get_lon_index(lon, lon_bounds, cross180)

            # extract lat/lon values
            lat_sub = lat[lat_ix]
            lon_sub = lon[lon_ix]

            xlon = track[:,2]

            z_sub = z[z_ix]

            if cross180:
                negative_lon = lon_sub < 0
                lon_sub[negative_lon] = lon_sub[negative_lon] + 360
                xlon[lix] = xlon[lix] + 360

            data = np.nan*np.ones((12, len(z_sub), len(lat_sub), len(lon_sub)))

        if verbose:
            sys.stdout.write('Extracting WOA data for {}\n'.format(mdates.num2date(mdates.datestr2num('2020-{:02d}-01'.format(i+1))).strftime('%b')))
        data[i,:,:,:] = nc.variables[var_name][:].data[0,z_ix,:,:][:,lat_ix,:][:,:,lon_ix]

    data[data > 9e36] = np.nan

    xtrack = track.copy()
    xtrack[:,2] = xlon
    woa_track = [z_sub, lat_sub, lon_sub]

    return xtrack, woa_track, data

def load_ncep_data(track, varname, local_path='./'):
    '''
    Function to load NCEP reanalysis data for comparison with autonomous
    float in-air data. Data to be interpolated along the provided
    track (t, lat, lon).

    INPUT:
              track: array with the columns (SDN, lat, lon)
              local_path: local directory where NCEP files are stored, assumes
                          current directory if no input

    OUTPUT:

    AUTHOR:   Christopher Gordon
              Fisheries and Oceans Canada
              chris.gordon@dfo-mpo.gc.ca

    LAST UPDATE: 04-05-2020

    CHANGE LOG:
    '''

    # make local_path a Path() object from a string, account for windows path
    local_path = Path(local_path)

    if varname == 'pres':
        base_file = 'pres.sfc.gauss'
        land_file = local_path / 'land' / 'land.sfc.gauss.nc'
    elif varname == 'rhum':
        base_file = 'rhum.sig995'
        land_file = local_path / 'land' / 'land.nc'
    else:
        raise ValueError('Invalid varname input')

    lnc = Dataset(land_file, 'r')

    # check if float track crosses -180/180 meridian
    cross180 = False
    if np.max(np.abs(np.diff(track[:,2]))) > 340:
        cross180 = True
        lix = track[:,2] < 0
        lon_bounds = (np.max(track[lix,2]), np.min(track[~lix,2]))
    else:
        lon_bounds = (np.min(track[:,2]), np.max(track[:,2]))
    lat_bounds = (np.min(track[:,1]), np.max(track[:,1]))

    sdn = track[:,0]
    yrs = (mdates.num2date(np.min(sdn)).year, mdates.num2date(np.max(sdn)).year)
    Nyear = yrs[1]-yrs[0]

    if Nyear == 0:
        Nyear = 1

    # counter index for going across years
    j = 0
    for n in range(Nyear+1):
        y = yrs[0] + n
        ncep_file = local_path / varname / '{}.{}.nc'.format(base_file, y)
        nc = Dataset(ncep_file, 'r')

        time = nc.variables['time'][:]
        time = time/24 + mdates.datestr2num('1800-01-01')

        if y == yrs[0]:
            lat = nc.variables['lat'][:]
            lon = nc.variables['lon'][:]
            lon[lon > 180] = lon[lon > 180] - 360
            lat_ix = util.get_lat_index(lat, lat_bounds)
            lon_ix = util.get_lon_index(lon, lon_bounds, cross180)

            # extract lat/lon values
            lat_sub = lat[lat_ix]
            lon_sub = lon[lon_ix]

            xlon = track[:,2]

            if cross180:
                negative_lon = lon_sub < 0
                lon_sub[negative_lon] = lon_sub[negative_lon] + 360
                xlon[lix] = xlon[lix] + 360

            extra_time_indices = 0
            for m in range(1,Nyear+1):
                ytest = yrs[0] + m
                ncep_file_test = local_path / varname / '{}.{}.nc'.format(base_file, ytest)
                nctest = Dataset(ncep_file_test, 'r')
                extra_time_indices += nctest.variables['time'].shape[0] - len(time)

            landmask = lnc.variables['land'][:][0,:,:][:,lon_ix][lat_ix,:].astype(bool)
            ncep_time = np.nan*np.ones((len(time)*(Nyear+1)+extra_time_indices))
            data = np.nan*np.ones((len(time)*(Nyear+1)+extra_time_indices, len(lat_sub), len(lon_sub)))

        vdata = nc.variables[varname][:]
        for i in range(len(time)):
            data_2d =  vdata[i,:,:][:,lon_ix][lat_ix,:]
            data_2d[landmask] = np.nan
            data[j,:,:] = data_2d
            ncep_time[j] = time[i]
            j += 1

        xtrack = track.copy()
        xtrack[:,2] = xlon
        ncep_track = [ncep_time, lat_sub, lon_sub]

    return xtrack, ncep_track, data
