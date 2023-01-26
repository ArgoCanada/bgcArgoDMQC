import sys
import warnings
from pathlib import Path

import numpy as np
import matplotlib.dates as mdates

from netCDF4 import Dataset

from .. import util

def load_woa_data(track, param, zlim=(0,1000), local_path='./', verbose=True):
    '''
    Function to load WOA18 climatological data for comparison with autonomous
    floats. Data to be interpolated along the provided track (t, lat, lon).

    Args:
        track: array with the columns (SDN, lat, lon)
        param: requested variable, valid inputs are
            - T: temperature
            - S: salinity
            - O2: dissolved oxygen
            - O2sat: oxygen percent saturation
            - NO3: nitrate
            - Si: silicate
            - PO4: phosphate
        zlim: depth bounds (upper, lower), default to (0, 1000)
        local_path: local directory where WOA files are stored, assumes
                    current directory if no input

    Returns:
        xtrack: same as track input, but adjusted lon if the track crosses the 180/-180 meridian
        woa_track: list with z, lat, and lon arrays of WOA data
        data: gridded array of the input variable (month, z, lat, lon)

    Acknowledgement: 
        This code is adapted from the SOCCOM SAGE_O2Argo matlab
        code, available via https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING,
        written by Tanya Maurer & Josh Plant
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
        lix = None
    lat_bounds = (np.nanmin(track[:,1]), np.nanmax(track[:,1]))

    # set up extraction files, variables
    woa_param, woa_ftype, woa_dir = util.decode_woa_var(param)
    var_name  = woa_param + '_an'

    base_woa_file = 'woa18_{}_{}'.format(woa_ftype, woa_param)

    # assign var names to avoid unbound warnings
    data, xlon, z_sub, lat_sub, lon_sub = None, None, None, None, None
    z_ix, lat_ix, lon_ix = None, None, None

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

    Args:
        track: array with the columns (SDN, lat, lon)
        local_path: local directory where NCEP files are stored, assumes current directory if no input

    Returns:

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
        lon_bounds = (np.nanmax(track[lix,2]), np.nanmin(track[~lix,2]))
    else:
        lon_bounds = (np.nanmin(track[:,2]), np.nanmax(track[:,2]))
        lix = None
    lat_bounds = (np.nanmin(track[:,1]), np.nanmax(track[:,1]))

    sdn = track[:,0]
    yrs = (mdates.num2date(np.nanmin(sdn)).year, mdates.num2date(np.nanmax(sdn)).year)
    Nyear = yrs[1]-yrs[0]

    if Nyear == 0 and yrs[0] != mdates.datetime.date.today().year:
        Nyear = 1

    # assign var names to avoid unbound warnings
    data, xlon, ncep_time, landmask, xtrack, ncep_track = None, None, None, None, None, None
    lat_sub, lon_sub, lat_ix, lon_ix = None, None, None, None

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
            data_2d = vdata[i,:,:][:,lon_ix][lat_ix,:]
            data_2d[landmask] = np.nan
            data[j,:,:] = data_2d
            ncep_time[j] = time[i]
            j += 1

        xtrack = track.copy()
        xtrack[:,2] = xlon
        ncep_track = [ncep_time, lat_sub, lon_sub]

    return xtrack, ncep_track, data

