from time import time
from pathlib import Path
import ftplib

import numpy as np

from .. import util
from .index import get_dac, URL_DIR_DICT, URL

def get_woa18(varname, local_path='./', ftype='netcdf', overwrite=False, nfiles=None):
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
    url = 'ftp-oceans.ncei.noaa.gov'
    param, dtype, ftpdir = util.decode_woa_var(varname)

    ftp = ftplib.FTP(url, 'anonymous', 'chris.gordon@dfo-mpo.gc.ca')
    ftp.cwd(f'pub/woa/WOA18/DATA/{ftpdir}/{ftype}/{dtype}/1.00/')

    local_path = local_path / ftpdir

    if not local_path.is_dir():
        local_path.mkdir()

    file_list = ftp.nlst()
    if nfiles  is not None and nfiles < len(file_list):
        file_list = file_list[:nfiles]

    for fn in file_list:
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
    else:
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

def get_argo(*args, local_path='./', url=URL, overwrite=False, summary_overwrite=True, ftype=None, mission='CB', mode='RD', nfiles=None):
    '''
    Function to download all data from a single float, or individual
    profiles

    Args:
        *args: list of files, floats, or directories
        local_path (optional, str or Path): local directory to save float data, defaults to current directory
        url (optional, str): url of the GDAC to connect to, defaults to ifremer or the value in the configuration file
        overwrite (optional, bool): whether to overwrite existing local files, default *False*
        ftype (optional, str): can be 'summary' if the user wishes to only download Sprof, traj, meta files, default is *None*
        mission (optional, str): Argo mission, can be 'B' for biogeochemical or 'C' for core, or 'CB' for both, default is 'CB'
        mode (optional, str): Download real-time ('R') or delayed-mode ('D') or all ('RD') data, default 'RD'
        
    Returns:
        ftplib.FTP object

    Author:   Christopher Gordon
              Fisheries and Oceans Canada
              chris.gordon@dfo-mpo.gc.ca

    Last update: 29-04-2020

    Change log:

    2020-10-27: add nfiles parameter to limit number of files downloaded,
    helpful to speed up testing/coverage runs
    '''

    if len(args) == 1:
        local_path = Path(local_path)

        ftp = ftplib.FTP(url)
        ftp.login()

        if type(args[0]) is not list:
            arglist = [args[0]]
        else:
            arglist = args[0]

        for init_a in arglist:
            # if its a float number, build ftp paths to floats
            if type(init_a) in [int, float, np.int32, np.int64, np.float32, np.float64]:
                base_path = URL_DIR_DICT[url]
                a = (Path(base_path) / 'dac' / get_dac(init_a) / str(init_a)).as_posix()
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
                    wmo_path = dac_path / wmo
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
                    # only download the file if 
                    # local_path (str or Path): local directory to save float data
                    # url (str): it doesn't already exist locally
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
                wmo_path = dac_path / wmo

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
                    if not wmo_file.exists() or overwrite or summary_overwrite:
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


                    if nfiles is not None and nfiles < len(files):
                        files = files[:nfiles]
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

    else:
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
                if not wmo_file.exists() or overwrite or summary_overwrite:
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

                if nfiles is not None and nfiles < len(files):
                    files = files[:nfiles]

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


