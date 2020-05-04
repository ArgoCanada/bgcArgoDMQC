#!/usr/bin/python

import sys
import warnings
from pathlib import Path

import numpy as np
import pylab as pl
from scipy.interpolate import interp1d

from netCDF4 import Dataset

from . import unit
from . import util

def apply_qc_adjustment():

    return None

def load_argo_data(local_path, wmo):
    # -------------------------------------------------------------------------
    # load_argo_data
    # -------------------------------------------------------------------------
    #
    # Function to load in all data from a single float, using BRtraj, meta,
    # and Sprof files
    #
    # INPUT:
    #           local_path: local path of float data
    #           wmo: float ID number
    #
    # OUTPUT:
    #           floatData: python dict() object with the following fields
    #               floatName: WMO number, from input
    #               floatType: Kind of float (APEX, ARVOR, etc.)
    #               N_LEVELS: Number of depth levels, Argo dimension N_LEVELS
    #               N_CYCLES: Number of profiles, Argo dimension N_PROF
    #               CYCLES: Array from 1 to N_CYCLES
    #               LATITUDE: Latitude (-90, 90) for each profile
    #               LONGITUDE: Longitude (-180, 180) for each profile
    #               SDN: Serial Date Number for each profile
    #               PRES: Pressure (dbar), compressed to vector (1D array)
    #               TEMP: Temperature (deg C)
    #               PSAL: Salinity (psu)
    #               DOXY: Dissolved Oxygen (micromole/kg)
    #               O2sat: Oxygen percent saturation (%)
    #               PPOX_DOXY: Oxygen partial pressure (atm) [if avail.]
    #               TRAJ_CYCLE: Cycle number for PPOX_DOXY [if avail.]
    #               inair: Boolean to indicate if in-air data exists
    #
    #               *** CYCLES, LATITUDE, LONGITUDE, AND SDN ALL ALSO HAVE ***
    #               ***     ANALOGOUS <VAR>_GRID FIELDS THAT MATCH THE     ***
    #               ***   DIMENSION OF PRES, TEMP, PSAL, DOXY, AND O2SAT   ***
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
    # 22-04-2020: updated so that pressure mask determines all variables - need
    # to add all quality flags to output
    #
    # 29-04-2020: switched file/path handling from os module to pathlib
    #
    # -------------------------------------------------------------------------

    # make local_path a Path() object from a string, account for windows path
    local_path = Path(local_path)

    # check that necessary files exist - can continue without BRtraj file but
    # need Sprof and meta files
    BRtraj = local_path / wmo / '{}_BRtraj.nc'.format(wmo)
    Sprof  = local_path / wmo / '{}_Sprof.nc'.format(wmo)
    meta   = local_path / wmo / '{}_meta.nc'.format(wmo)

    # check if BRtraj is there, flag for moving forward if not
    BRtraj_flag = True
    if not BRtraj.exists():
        BRtraj_flag = False
        sys.stdout.write('Continuing without BRtraj file\n')

    # Sprof and meta are required, so raise error if they are not there
    if not Sprof.exists():
        raise FileNotFoundError('No such Sprof file: {}'.format(Sprof))
    if not meta.exists():
        raise FileNotFoundError('No such meta file: {}'.format(meta))

    # only load if file was found
    if BRtraj_flag:
        BRtraj_nc = Dataset(BRtraj)
    else:
        BRtraj_nc = None
    # load synthetic and meta profiles
    Sprof_nc = Dataset(Sprof)
    meta_nc  = Dataset(meta)

    # number of profile cycles
    M = Sprof_nc.dimensions['N_LEVELS'].size
    N = Sprof_nc.dimensions['N_PROF'].size
    # beginning of output dict with basic info, following variables in SAGEO2
    floatData = dict(floatName=wmo, N_CYCLES=N, N_LEVELS=M, CYCLES=np.arange(1,N+1))

    ftype = ''
    for let in meta_nc.variables['PLATFORM_TYPE'][:].compressed():
        ftype = ftype + let.decode('UTF-8')

    floatData['floatType'] = ftype

    pres = Sprof_nc.variables['PRES'][:]

    t = Sprof_nc.variables['JULD'][:].compressed() + pl.datestr2num('1950-01-01')
    t_grid = np.ma.masked_array(np.tile(t,(M,1)).T, mask=pres.mask)

    lat = Sprof_nc.variables['LATITUDE'][:].compressed()
    lon = Sprof_nc.variables['LONGITUDE'][:].compressed()

    lat_grid = np.ma.masked_array(np.tile(lat,(M,1)).T, mask=pres.mask).compressed()
    lon_grid = np.ma.masked_array(np.tile(lon,(M,1)).T, mask=pres.mask).compressed()

    cycle_grid = np.ma.masked_array(np.tile(floatData['CYCLES'],(M,1)).T, mask=pres.mask)

    # use the pressure mask for all variables to ensure dimensions match
    floatData['PRES'] = pres.compressed()
    floatData['TEMP'] = np.ma.masked_array(Sprof_nc.variables['TEMP'][:].data, mask=pres.mask).compressed()
    floatData['PSAL'] = np.ma.masked_array(Sprof_nc.variables['PSAL'][:].data, mask=pres.mask).compressed()
    floatData['DOXY'] = np.ma.masked_array(Sprof_nc.variables['DOXY'][:].data, mask=pres.mask).compressed()

    floatData['SDN'] = t
    floatData['SDN_GRID']  = t_grid.compressed()
    floatData['CYCLE_GRID'] = cycle_grid.compressed()

    floatData['LATITUDE'] = lat
    floatData['LONGITUDE'] = lon
    floatData['LATITUDE_GRID'] = lat_grid
    floatData['LONGITUDE_GRID'] = lon_grid

    floatData['O2Sat'] = 100*floatData['DOXY']/unit.oxy_sol(floatData['PSAL'], floatData['TEMP'], unit='micromole/kg')

    if BRtraj_flag:
        ppox_doxy = BRtraj_nc.variables['PPOX_DOXY'][:]
        floatData['pO2'] = ppox_doxy.compressed()
        floatData['TRAJ_CYCLE'] = np.ma.masked_array(BRtraj_nc.variables['CYCLE_NUMBER'][:].data, mask=ppox_doxy.mask).compressed()
        floatData['inair'] = True
    else:
        floatData['inair'] = False

    return floatData

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

def interp_woa_data(track, woa_track, data):
    # -------------------------------------------------------------------------
    # load_woa_data
    # -------------------------------------------------------------------------
    #
    # Function to interpolate WOA18 climatological data along the provided 
    # track (t, lat, lon).
    #
    # INPUT:
    #           track: array with the columns (SDN, lat, lon)
    #           woa_track: tuple with WOA time, lat, lon arrays
    #           data: output array from load_woa_data()
    #
    # OUTPUT:
    #           woa_interp: 2D array of requested WOA parameter (depth x time)
    #
    # AUTHOR:   Christopher Gordon
    #           Fisheries and Oceans Canada
    #           chris.gordon@dfo-mpo.gc.ca
    #
    # ACKNOWLEDGEMENT: this code is adapted from the SOCCOM SAGE_O2Argo matlab
    # code, available via https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING,
    # written by Tanya Maurer & Josh Plant
    #
    # LAST UPDATE: 23-04-2020
    #
    # CHANGE LOG:
    #
    # -------------------------------------------------------------------------

    # get relevant WOA information
    z, lat, lon = woa_track

    # put float and WOA on common time axis/year
    M = z.shape[0]
    N = track.shape[0]
    yrday = np.array([dn - pl.datestr2num('{}-01-01'.format(pl.num2date(dn).year)) for dn in track[:,0]])
    woa_yrday = np.array([pl.datestr2num('2020-{:02d}-15'.format(i+1)) - pl.datestr2num('2020-01-01') for i in range(12)])

    # array for output
    woa_interp = np.nan*np.ones((M,N))

    xwt = [[],[],[]]

    for i in range(N):
        # ---------------------------------------------------------------------
        # Get bounding indices, interp weights, and subset before interpolation
        # ---------------------------------------------------------------------
        t1 = np.where(woa_yrday < yrday[i])[0]
        if t1.shape[0] == 0: # early Jan
            t1 = 0
            t2 = 11
            wt = (yrday[i] + 15) / 30
        elif t1[-1] == 11: # late Dec
            t1 = t1[-1]
            t2 = 0
            wt = (365 - yrday[i] + 15) / 30
        else:
            t1 = t1[-1]
            t2 = t1 + 1
            dx1 = woa_yrday[t2] - woa_yrday[t1]
            dx2 = yrday[i] - woa_yrday[t1]
            wt = (dx1 - dx2) / dx1

        xwt[0].append(wt)

        lat_ix1 = np.where(lat < track[i,1])[0][-1]
        lat_ix2 = lat_ix1 + 1
        dx1 = lat[lat_ix2] - lat[lat_ix1]
        dx2 = track[i,1] - lat[lat_ix1]
        lat_wt = (dx1 - dx2) / dx1

        xwt[1].append(lat_wt)

        lon_ix1 = np.where(lon < track[i,2])[0][-1]
        lon_ix2 = lon_ix1 + 1
        dx1 = lon[lon_ix2] - lon[lon_ix1]
        dx2 = track[i,2] - lon[lon_ix1]
        lon_wt = (dx1 - dx2) / dx1

        xwt[2].append(lon_wt)

        # ---------------------------------------------------------------------
        # Interpolation part 1 - get bounding profiles, interp over time
        # ---------------------------------------------------------------------
        D3 = wt*data[t1,:,lat_ix1:lat_ix2+1,lon_ix1:lon_ix2+1] + (1-wt)*data[t2,:,lat_ix1:lat_ix2+1,lon_ix1:lon_ix2+1]

        # ---------------------------------------------------------------------
        # Interpolation part 2 - interp over lat + lon
        # ---------------------------------------------------------------------
        if np.isnan(D3.flatten()).any():
            sys.stdout.write('Bounding climatological profile(s) missing data')
            sys.stdout.write(' - taking simple average of available data.\n')
            woa_interp[:,i] = np.nanmean(np.nanmean(D3, axis=2), axis=1)
        else:
            D2 = lat_wt*D3[:,0,:] + (1 - lat_wt)*D3[:,1,:]
            woa_interp[:,i] = lon_wt*D2[:,0] + (1 - lon_wt)*D2[:,1]

    return woa_interp, xwt, yrday

def woa_to_float_track(track, param, zlim=(0,1000), local_path='./'):
    # -------------------------------------------------------------------------
    # woa_to_float_track
    # -------------------------------------------------------------------------
    #
    # Function to load WOA18 climatological data for comparison with autonomous
    # floats. Data to be interpolated along the provided track (t, lat, lon).
    # Combines function load_woa_data() and interp_woa_data() for convenience,
    # see documentation for those funcions for more detail.
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
    #           z: WOA depth array
    #           woa_interp: 2D array of requested WOA parameter (depth x time)
    #
    # AUTHOR:   Christopher Gordon
    #           Fisheries and Oceans Canada
    #           chris.gordon@dfo-mpo.gc.ca
    #
    # LAST UPDATE: 23-04-2020
    #
    # CHANGE LOG:
    #
    # -------------------------------------------------------------------------

    xtrack, woa_track, woa_data = load_woa_data(track, param, zlim=zlim, local_path=local_path)
    woa_interp, wt, yrday = interp_woa_data(xtrack, woa_track, woa_data)
    z = woa_track[0]

    return z, woa_interp

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

def interp_ncep_data(track, ncep_track, data):

    # extract ncep variables
    ncep_time, lat_sub, lon_sub = ncep_track

    # get float lats and lons, round out times at endpoints
    xt   = np.ones((track.shape[0]+2))
    xlat = np.ones((track.shape[0]+2))
    xlon = np.ones((track.shape[0]+2))

    xt[1:-1]   = track[:,0]
    xlat[1:-1] = track[:,1]
    xlon[1:-1] = track[:,2]

    xt[0]   = np.floor(xt[1])
    xt[-1]  = np.ceil(xt[-2])
    xlat[0] = xlat[1]
    xlat[-1] = xlat[-2]
    xlon[0] = xlon[1]
    xlon[-1] = xlon[-2]

    # interpolate float to ncep time grid
    f = interp1d(xt, xlat, bounds_error=False)
    ilat = f(ncep_time)
    f = interp1d(xt, xlon, bounds_error=False)
    ilon = f(ncep_time)

    # remove NaN values
    ix = ~np.isnan(ilat)
    ilat = ilat[ix]
    ilon = ilon[ix]
    time_sub = ncep_time[ix]

    # loop through NCEP time and get data at each location
    ncep_data = np.nan*np.ones((time_sub.shape[0],))
    for i in range(time_sub.shape[0]):
        # temporal weighted average
        time_pt = np.where(ncep_time >= time_sub[i])[0][0]
        if time_pt == 0:
            tmp = data[time_pt,:,:]
        else:
            wt = (ncep_time[time_pt] - time_sub[i]) / (ncep_time[time_pt] - ncep_time[time_pt-1])
            tmp = wt*data[time_pt-1,:,:] + (1-wt)*data[time_pt,:,:]
        
        # latitudinal weighted average
        lat_pt = np.where(lat_sub >= ilat[i])[0][-1]
        wt = (lat_sub[lat_pt] - ilat[i]) / (lat_sub[lat_pt] - lat_sub[lat_pt+1])
        tmp = wt*tmp[lat_pt+1,:] + (1-wt)*tmp[lat_pt,:]

        # longitudinal weighted average - final pt. 
        lon_pt = np.where(lon_sub >= ilon[i])[0][0]
        wt = (lon_sub[lon_pt] - ilon[i]) / (lon_sub[lon_pt] - lon_sub[lon_pt-1])
        tmp = wt*tmp[lon_pt-1] + (1-wt)*tmp[lon_pt]

        ncep_data[i] = tmp

    # interpolate back to float time
    f = interp1d(time_sub, ncep_data, bounds_error=False)
    ncep_interp = f(track[:,0])

    return ncep_interp

def ncep_to_float_track():

    return None


def calc_gain(data, ref, inair=True, zlim=25.):
    # -------------------------------------------------------------------------
    # calc_gain
    # -------------------------------------------------------------------------
    #
    # Calculate the gain for each profile by comparing float oxygen data to a
    # reference data set, either NCEP for in-air or WOA surface data if in-air
    # comparison is not available.
    #
    # INPUT:
    #           data: float data dict object, output from load_argo_data()
    #           ref: reference data set, either NCEP pO2 or WOA O2sat
    #           inair: boolean flag to indicate if comparison to NCEP in-air
    #                  data or WOA surface data should be done, default to
    #                  in-air, but function also performs check
    #           zlim: lower limit to define as 'surface' and take mean within,
    #                 default value 25 dbar
    #
    # OUTPUT:
    #           g: vector of gains
    #           surf_data: array of float surface stats (cycle, N, mean, std)
    #
    # AUTHOR:   Christopher Gordon
    #           Fisheries and Oceans Canada
    #           chris.gordon@dfo-mpo.gc.ca
    #
    # LAST UPDATE: 23-04-2020
    #
    # CHANGE LOG:
    #
    # -------------------------------------------------------------------------

    # check which reference data to use
    if inair and 'pO2' not in data.keys():
        raise ValueError('Flag ''inair'' set to True but partial pressure data not available')

    if inair:
        sys.stdout.write('NCEP/in-air functionality not built yet, returning array of ones for sake of testing\n')
        g = np.nan*np.ones((ref.shape[0],))

        surf_ix = data['PRES'] <= zlim
        surf_ppox = data['O2Sat'][surf_ix]
        cycle = data['CYCLE_GRID'][surf_ix]

        surf_data = np.nan*np.ones((ref.shape[0],4))
        for i,c in enumerate(np.unique(cycle)):
            subset_ppox = surf_ppox[cycle == c]
            surf_data[i,0] = c
            surf_data[i,1] = np.sum(~np.isnan(subset_ppox))
            surf_data[i,2] = np.nanmean(subset_ppox)
            surf_data[i,3] = np.nanstd(subset_ppox)

            g[i] = ref[i]/surf_data[i,2]
    else:
        sys.stdout.write('Calculating gains using WOA surface data and float O2 percent saturation')
        surf_ix = data['PRES'] <= zlim
        surf_o2sat = data['O2Sat'][surf_ix]
        cycle = data['CYCLE_GRID'][surf_ix]

        z_woa = ref['z']
        woa_data = ref['data']

        woa_index = np.where(z_woa <= zlim)[0]
        woa_surf = np.nanmean(woa_data[woa_index,:],axis=0)
        woa_surf = woa_data[0,:]

        surf_data = np.nan*np.ones((woa_surf.shape[0],4))
        g = np.nan*np.ones((woa_surf.shape[0],))
        for i,c in enumerate(np.unique(cycle)):
            ref_o2sat = woa_surf[i]
            subset_o2sat = surf_o2sat[cycle == c]
            surf_data[i,0] = c
            surf_data[i,1] = np.sum(~np.isnan(subset_o2sat))
            surf_data[i,2] = np.nanmean(subset_o2sat)
            surf_data[i,3] = np.nanstd(subset_o2sat)

            g[i] = ref_o2sat/surf_data[i,2]
        

    return g, surf_data

def aic(data,resid):
    # -------------------------------------------------------------------------
    # aic
    # -------------------------------------------------------------------------
    #
    # function to calculate the Akiake Information Criteria (AIC) as a metric
    # for assessing the appropriate number of breakpoints in the calculation of
    # drifts in O2 gains.
    #
    # INPUT:
    #
    # OUTPUT:
    #
    # AUTHOR:   Christopher Gordon
    #           Fisheries and Oceans Canada
    #           chris.gordon@dfo-mpo.gc.ca
    #
    # ACKNOWLEDGEMENT: this code is adapted from the SOCCOM SAGE_O2Argo matlab
    # code, available via https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING,
    # written by Tanya Maurer & Josh Plant
    #
    # LAST UPDATE: 20-04-2020
    #
    # CHANGE LOG:
    #
    # -------------------------------------------------------------------------

    # calculate AIC
    SSE = np.sum(resid**2) # sum square errors
    n = resid.shape[0]
    m = data.shape[0] - 1 # do not include first cycle
    K = 2*m + 2

    # valid data parameters? see Jones & Day (1995)
    is_valid = n/4 - 1
    if m > is_valid:
        aic_value = np.nan
        sys.stdout.write('n >> K, cannot caclculate AIC, setting AIC = NaN')
    else:
        # formula ref. Jones & Day (1995), Owens & Wong (2009)
        aic_value = np.log(SSE/n) + (n+K)/(n-K-2)

    return aic_value


def bic(data,resid):
    # -------------------------------------------------------------------------
    # bic
    # -------------------------------------------------------------------------
    #
    # function to calculate the Bayesian Information Criteria (BIC) as a metric
    # for assessing the appropriate number of breakpoints in the calculation of
    # drifts in O2 gains.
    #
    # INPUT:
    #
    # OUTPUT:
    #
    # AUTHOR:   Christopher Gordon
    #           Fisheries and Oceans Canada
    #           chris.gordon@dfo-mpo.gc.ca
    #
    # ACKNOWLEDGEMENT: this code is adapted from the SOCCOM SAGE_O2Argo matlab
    # code, available via https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING,
    # written by Tanya Maurer & Josh Plant
    #
    # LAST UPDATE: 20-04-2020
    #
    # CHANGE LOG:
    #
    # -------------------------------------------------------------------------

    # calculate BIC
    errorlim = 0 # cap on residuals, useful for noisy pH and nitrate data
    SSE = np.sum(resid**2) # sum square errors
    n = resid.shape[0]
    m = data.shape[0] - 1 # do not include first cycle
    K = 2*m + 2

    # valid data parameters? see Jones & Day (1995)
    is_valid = n/4 - 1
    if m > is_valid:
        bic_value = np.nan
        sys.stdout.write('n >> K, cannot caclculate BIC, setting BIC = NaN')
    else:
        bic_value = np.log(1/(n*SSE) + errorlim**2) + K*np.log(n)/n

    return bic_value
