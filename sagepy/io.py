#!/usr/bin/python

import sys

from pathlib import Path
import ftplib
import gzip

from . import util

def get_woa18(varname, local_path='./', ftype='netcdf', overwrite=False):
    # -------------------------------------------------------------------------
    # get_woa18
    # -------------------------------------------------------------------------
    #
    # Function to download WOA data for a given variable
    #
    # INPUT:
    #           varname: woa18 variable to download data for - one of:
    #               T: temperature
    #               S: salinity
    #               O2: dissolved oxygen
    #               O2sat: oxygen percent saturation
    #               NO3: nitrate
    #               Si: silicate
    #               PO4: phosphate
    #           local_path: path to save files to (inside variable folder, ex.
    #               if local_path = '../data/woa', and varname = 'T', files
    #               will be saved in ../data/woa/temperature/*.nc. Defaults
    #               to current directory
    #           ftype: file format to download, defaults to netcdf (.nc):
    #               ascii: .dat
    #               csv: .csv
    #               netcdf: .nc
    #           overwrite: boolean flag, if False, does not re-download
    #               existing files, if true, will download no matter what, 
    #               defaults to False
    #
    # OUTPUT:
    #           ftp: the ftplib server object
    #
    # AUTHOR:   Christopher Gordon
    #           Fisheries and Oceans Canada
    #           chris.gordon@dfo-mpo.gc.ca
    #
    # LAST UPDATE: 29-04-2020
    #
    # CHANGE LOG:
    #
    # -------------------------------------------------------------------------

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
        if not local_file.exists() | overwrite:
            print(local_file)
            # open the local file
            lf = open(local_file, 'wb')
            # retrieve the file on FTP server,
            ftp.retrbinary('RETR ' + fn, lf.write)


    return ftp

def get_ncep(varname, local_path='./', overwrite=False):
    # -------------------------------------------------------------------------
    # get_woa18
    # -------------------------------------------------------------------------
    #
    # Function to download NCEP reanalysis gaussian gridded surface air 
    # pressure data 
    #
    # INPUT:
    #           varname: 'pres' (pressure) or 'rhum' (relative humidity)
    #               or 'land' (to get land mask)
    #           local_path: path to save files to, defaults
    #               to current directory
    #           overwrite: boolean flag, if False, does not re-download
    #               existing files, if true, will download no matter what, 
    #               defaults to False
    #
    # OUTPUT:
    #           ftp: the ftplib server object
    #
    # AUTHOR:   Christopher Gordon
    #           Fisheries and Oceans Canada
    #           chris.gordon@dfo-mpo.gc.ca
    #
    # LAST UPDATE: 29-04-2020
    #
    # CHANGE LOG:
    #
    # -------------------------------------------------------------------------

    local_path = Path(local_path)
    url = 'ftp.cdc.noaa.gov'

    ftp = ftplib.FTP(url, 'anonymous', 'chris.gordon@dfo-mpo.gc.ca')

    if varname == 'pres':
        ftp.cwd('Datasets/ncep.reanalysis2/gaussian_grid/')

        local_path = local_path / varname
        if not local_path.is_dir():
            local_path.mkdir()

        for yr in range(2010, 2021):
            fn = 'pres.sfc.gauss.{}.nc'.format(yr)
            local_file = local_path / fn

            if not local_file.exists() | overwrite:
                print(local_file)
                # open the local file
                lf = open(local_file, 'wb')
                # retrieve the file on FTP server,
                ftp.retrbinary('RETR ' + fn, lf.write)
    
    elif varname == 'rhum':

        ftp.cwd('Datasets/ncep.reanalysis/surface/')

        local_path = local_path / varname
        if not local_path.is_dir():
            local_path.mkdir()

        for yr in range(2010, 2021):
            fn = 'rhum.sig995.{}.nc'.format(yr)
            local_file = local_path / fn
            if not local_file.exists() | overwrite:
                print(local_file)
                # open the local file
                lf = open(local_file, 'wb')
                # retrieve the file on FTP server,
                ftp.retrbinary('RETR ' + fn, lf.write)

    elif varname == 'land':

        local_path = local_path / varname
        if not local_path.is_dir():
            local_path.mkdir()

        ftp.cwd('Datasets/ncep.reanalysis2/gaussian_grid/')
        fn = 'land.sfc.gauss.nc'
        local_file = local_path / fn
        if not local_file.exists() | overwrite:
            lf = open(local_file, 'wb')
            ftp.retrbinary('RETR ' + fn, lf.write)

        ftp.cwd('../../ncep.reanalysis/surface/')
        fn = 'land.nc'
        local_file = local_path / fn
        if not local_file.exists() | overwrite:
            lf = open(local_file, 'wb')
            ftp.retrbinary('RETR ' + fn, lf.write)
    
    else:
        raise ValueError('Invalid varname input')

    return ftp

def get_argo(*args, local_path='./', url='ftp.ifremer.fr', overwrite=False):
    # -------------------------------------------------------------------------
    # get_argo
    # -------------------------------------------------------------------------
    #
    # Function to download all data from a single float, or individual 
    # profiles
    #
    # INPUT:
    #           Inputs may vary depending on desired performance. Multiple
    #           arguments may be provided to download all files from a certain
    #           float or argo defined geographical area. A single path to a 
    #           file may be provided to download that file. A list of files 
    #           may be provided as well.
    #
    #           overwrite: boolean flag, if False, does not re-download
    #               existing files, if true, will download no matter what, 
    #               defaults to False
    #
    # OUTPUT:
    #
    # AUTHOR:   Christopher Gordon
    #           Fisheries and Oceans Canada
    #           chris.gordon@dfo-mpo.gc.ca
    #
    # LAST UPDATE: 29-04-2020
    #
    # CHANGE LOG:
    #
    # -------------------------------------------------------------------------

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
            if not wmo_file.exists() | overwrite:
                print(wmo_file)
                # open the local file
                lf = open(wmo_file, 'wb')
                # retrieve the file on FTP server,
                ftp.retrbinary('RETR ' + fn, lf.write)

        # ------------------------ INDIVIDUAL PROFILE FILES -------------------
        # repeat as above
        if 'profiles' in ftp.nlst():
            ftp.cwd('profiles')
            files = ftp.nlst('*.nc')

            profile_path = wmo_path / 'profiles'
            if not profile_path.exists():
                profile_path.mkdir()

            for fn in files:
                profile_file = profile_path / fn
                if not profile_file.exists() | overwrite:
                    print(profile_file)
                    lf = open(profile_file, 'wb')
                    ftp.retrbinary('RETR ' + fn, lf.write)

            # back to parent directory
            ftp.cwd('../../')
        else:
            ftp.cwd('../')


    return ftp
