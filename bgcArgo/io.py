#!/usr/bin/python

import sys

from pathlib import Path
import ftplib
import gzip

import numpy as np
import pylab as pl

from netCDF4 import Dataset

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

def load_woa_data(track, param, zlim=(0,1000), local_path='./'):
    # -------------------------------------------------------------------------
    # load_woa_data
    # -------------------------------------------------------------------------
    #
    # Function to load WOA18 climatological data for comparison with autonomous
    # floats. Data to be interpolated along the provided track (t, lat, lon).
    #
    # INPUT:
    #           track: array with the columns (SDN, lat, lon)
    #           param: requested variable, valid inputs are
    #               T: temperature
    #               S: salinity
    #               O2: dissolved oxygen
    #               O2sat: oxygen percent saturation
    #               NO3: nitrate
    #               Si: silicate
    #               PO4: phosphate
    #           zlim: depth bounds (upper, lower), default to (0, 1000)
    #           local_path: local directory where WOA files are stored, assumes
    #                       current directory if no input
    #
    # OUTPUT:
    #           xtrack: same as track input, but adjusted lon if the track
    #                   crosses the 180/-180 meridian
    #           woa_track: list with z, lat, and lon arrays of WOA data
    #           data: gridded array of the input variable (month, z, lat, lon)
    #
    # AUTHOR:   Christopher Gordon
    #           Fisheries and Oceans Canada
    #           chris.gordon@dfo-mpo.gc.ca
    #
    # ACKNOWLEDGEMENT: this code is adapted from the SOCCOM SAGE_O2Argo matlab
    # code, available via https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING,
    # written by Tanya Maurer & Josh Plant
    #
    # LAST UPDATE: 29-04-2020
    #
    # CHANGE LOG:
    #
    # 23-04-2020: changed zlim to optional input argument
    #
    # 29-04-2020: switched file/path handling from os module to pathlib
    #
    # -------------------------------------------------------------------------

    # make local_path a Path() object from a string, account for windows path
    local_path = Path(local_path)

    # check if float track crosses -180/180 meridian
    cross180 = False
    if np.max(np.abs(np.diff(track[:,2]))) > 340:
        cross180 = True
        lix = track[:,2] < 0
        lon_bounds = (np.max(track[lix,2]), np.min(track[~lix,2]))
    else:
        lon_bounds = (np.min(track[:,2]), np.max(track[:,2]))
    lat_bounds = (np.min(track[:,1]), np.max(track[:,1]))

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

        sys.stdout.write('Extracting WOA data for {}\n'.format(pl.num2date(pl.datestr2num('2020-{:02d}-01'.format(i+1))).strftime('%b')))
        data[i,:,:,:] = nc.variables[var_name][:].data[0,z_ix,:,:][:,lat_ix,:][:,:,lon_ix]

    data[data > 9e36] = np.nan

    xtrack = track.copy()
    xtrack[:,2] = xlon
    woa_track = [z_sub, lat_sub, lon_sub]

    return xtrack, woa_track, data

def load_ncep_data(track, varname, local_path='./'):
    # -------------------------------------------------------------------------
    # load_ncep_data
    # -------------------------------------------------------------------------
    #
    # Function to load NCEP reanalysis data for comparison with autonomous
    # float in-air data. Data to be interpolated along the provided 
    # track (t, lat, lon).
    #
    # INPUT:
    #           track: array with the columns (SDN, lat, lon)
    #           local_path: local directory where NCEP files are stored, assumes
    #                       current directory if no input
    #
    # OUTPUT:
    #
    # AUTHOR:   Christopher Gordon
    #           Fisheries and Oceans Canada
    #           chris.gordon@dfo-mpo.gc.ca
    #
    # LAST UPDATE: 04-05-2020
    #
    # CHANGE LOG:
    #
    # -------------------------------------------------------------------------

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
    yrs = (pl.num2date(np.min(sdn)).year, pl.num2date(np.max(sdn)).year)
    Nyear = yrs[1]-yrs[0]
    
    # counter index for going across years
    j = 0
    for y in range(yrs[0], yrs[1]):
        ncep_file = local_path / varname / '{}.{}.nc'.format(base_file, y)
        nc = Dataset(ncep_file, 'r')

        time = nc.variables['time'][:]
        time = time/24 + pl.datestr2num('1800-01-01')

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
            
            landmask = lnc.variables['land'][:][0,:,:][:,lon_ix][lat_ix,:].astype(bool)
            ncep_time = np.nan*np.ones((len(time)*Nyear))
            data = np.nan*np.ones((len(time)*Nyear, len(lat_sub), len(lon_sub)))
        
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