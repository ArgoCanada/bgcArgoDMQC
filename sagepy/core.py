#!/usr/bin/python

import sys
import warnings
from pathlib import Path

import numpy as np
import pylab as pl

from netCDF4 import Dataset

from . import io
from . import interp
from . import unit
from . import util

def apply_qc_adjustment():

    return None

def argo(local_path, wmo):
    # -------------------------------------------------------------------------
    # argo
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

    if type(wmo) is not str:
        wmo = str(wmo)

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
    elif BRtraj.exists():
        BRtraj_nc = Dataset(BRtraj, 'r')
        if 'PPOX_DOXY' not in BRtraj_nc.variables.keys():
            BRtraj_flag = False
            sys.stdout.write('BRtraj file exists, but no in-air data exists, continuing without using BRtraj file\n')
    else:
        BRtraj_nc = None

    # Sprof and meta are required, so raise error if they are not there
    if not Sprof.exists():
        raise FileNotFoundError('No such Sprof file: {}'.format(Sprof))
    if not meta.exists():
        raise FileNotFoundError('No such meta file: {}'.format(meta))

    # load synthetic and meta profiles
    Sprof_nc = Dataset(Sprof, 'r')
    meta_nc  = Dataset(meta, 'r')

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
        floatData['PPOX_DOXY'] = ppox_doxy.compressed()
        floatData['TEMP_DOXY'] = np.ma.masked_array(BRtraj_nc.variables['TEMP_DOXY'][:].data, mask=ppox_doxy.mask).compressed()
        floatData['TRAJ_CYCLE'] = np.ma.masked_array(BRtraj_nc.variables['CYCLE_NUMBER'][:].data, mask=ppox_doxy.mask).compressed()
        floatData['inair'] = True
    else:
        floatData['inair'] = False

    return floatData
    
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

    xtrack, woa_track, woa_data = io.load_woa_data(track, param, zlim=zlim, local_path=local_path)
    woa_interp, wt, yrday = interp.interp_woa_data(xtrack, woa_track, woa_data)
    z = woa_track[0]

    return z, woa_interp

def ncep_to_float_track(varname, track, local_path='./'):

    xtrack, ncep_track, data = io.load_ncep_data(track, varname, local_path=local_path)
    ncep_interp = interp.interp_ncep_data(xtrack, ncep_track, data)

    return ncep_interp


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
    #               data or WOA surface data should be done, default to
    #               in-air, but function also performs check
    #           zlim: lower limit to define as 'surface' and take mean within,
    #                 default value 25 dbar, for use only when inair is False
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
    if inair and 'PPOX_DOXY' not in data.keys():
        raise ValueError('Flag ''inair'' set to True but partial pressure data not available')

    if inair:
        g = np.nan*np.ones((ref.shape[0],))

        # float partial pressure measurements at each cycle
        ppox  = data['PPOX_DOXY']
        cycle = data['CYCLES']

        mean_float_data = np.nan*np.ones((ref.shape[0],4))
        for i,c in enumerate(cycle):
            subset_ppox = ppox[data['TRAJ_CYCLE'] == c]
            mean_float_data[i,0] = c
            mean_float_data[i,1] = np.sum(~np.isnan(subset_ppox))
            mean_float_data[i,2] = np.nanmean(subset_ppox)
            mean_float_data[i,3] = np.nanstd(subset_ppox)

            g[i] = ref[i]/mean_float_data[i,2]

    else:
        sys.stdout.write('Calculating gains using WOA surface data and float O2 percent saturation')
        surf_ix = data['PRES'] <= zlim
        surf_o2sat = data['O2Sat'][surf_ix]
        cycle = data['CYCLE_GRID'][surf_ix]

        z_woa = ref['z']
        woa_data = ref['WOA']

        woa_index = np.where(z_woa <= zlim)[0]
        woa_surf = np.nanmean(woa_data[woa_index,:],axis=0)
        woa_surf = woa_data[0,:]

        mean_float_data = np.nan*np.ones((woa_surf.shape[0],4))
        g = np.nan*np.ones((woa_surf.shape[0],))
        for i,c in enumerate(np.unique(cycle)):
            ref_o2sat = woa_surf[i]
            subset_o2sat = surf_o2sat[cycle == c]
            mean_float_data[i,0] = c
            mean_float_data[i,1] = np.sum(~np.isnan(subset_o2sat))
            mean_float_data[i,2] = np.nanmean(subset_o2sat)
            mean_float_data[i,3] = np.nanstd(subset_o2sat)

            g[i] = ref_o2sat/mean_float_data[i,2]
        

    return g, mean_float_data

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

def get_var_by(v1, v2, float_data):
    
    index = np.unique(float_data[v2])
    out_array = np.nan*np.ones((len(index)))
    for i,v in enumerate(index):
        out_array[i] = np.nanmean(float_data[v1][float_data[v2] == v])
    
    return out_array
