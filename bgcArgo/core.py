#!/usr/bin/python

import sys
import warnings
from pathlib import Path
import fnmatch

import numpy as np
import pylab as pl
from scipy.interpolate import interp1d, RectBivariateSpline

import matplotlib.pyplot as plt

import seawater as sw

from netCDF4 import Dataset

from . import io
from . import interp
from . import unit
from . import util
from . import fplt
from . import diagnostic

# ----------------------------------------------------------------------------
# LOCAL MACHINE SETUP
# ----------------------------------------------------------------------------

ARGO_PATH = './'
WOA_PATH  = None
NCEP_PATH = None

__bgcindex__ = io.read_index()

def set_dirs(argo_path=ARGO_PATH, woa_path=WOA_PATH, ncep_path=NCEP_PATH):

    global ARGO_PATH
    ARGO_PATH = argo_path
    global WOA_PATH
    WOA_PATH  = woa_path
    global NCEP_PATH
    NCEP_PATH = ncep_path

# ----------------------------------------------------------------------------
# FLOAT CLASS
# ----------------------------------------------------------------------------

def get_index(index='bgc'):
    if index == 'bgc':
        return __bgcindex__
    elif index == 'global':
        return io.read_index(mission='C')
    elif index == 'synthetic':
        return io.read_index(mission='S')

class sprof:

    set_dirs = set_dirs

    def __init__(self, wmo):
        self.__floatdict__ = load_argo(ARGO_PATH, wmo, grid=True)

        # local path info
        self.argo_path = ARGO_PATH
        self.woa_path  = WOA_PATH
        self.ncep_path = NCEP_PATH

        # metadata and dimension variables
        self.floatName = self.__floatdict__['floatName']
        # self.floatType = self.__floatdict__['floatType']
        self.N_CYCLES  = self.__floatdict__['N_CYCLES']
        self.N_LEVELS  = self.__floatdict__['N_LEVELS']
        self.CYCLES    = self.__floatdict__['CYCLES']

        # time and location data
        self.SDN       = self.__floatdict__['SDN']
        self.LATITUDE  = self.__floatdict__['LATITUDE']
        self.LONGITUDE = self.__floatdict__['LONGITUDE']

        # core variables
        self.PRES = self.__floatdict__['PRES']
        self.TEMP = self.__floatdict__['TEMP']
        self.PSAL = self.__floatdict__['PSAL']
        # potential density
        self.PDEN = sw.pden(self.PSAL, self.TEMP, self.PRES) - 1000

        # bgc variables - not necessarily all there so check if the fields exist
        if 'DOXY' in self.__floatdict__.keys():
            self.DOXY = self.__floatdict__['DOXY']
        if 'CHLA' in self.__floatdict__.keys():
            self.CHLA = self.__floatdict__['CHLA']
        if 'BBP700' in self.__floatdict__.keys():
            self.BBP700 = self.__floatdict__['BBP700']
        if 'CDOM' in self.__floatdict__.keys():
            self.CDOM = self.__floatdict__['CDOM']
        if 'DOXY_ADJUSTED' in self.__floatdict__.keys():
            self.DOXY_ADJUSTED = self.__floatdict__['DOXY_ADJUSTED']
        if 'CHLA_ADJUSTED' in self.__floatdict__.keys():
            self.CHLA_ADJUSTED = self.__floatdict__['CHLA_ADJUSTED']
        if 'BBP700_ADJUSTED' in self.__floatdict__.keys():
            self.BBP700_ADJUSTED = self.__floatdict__['BBP700_ADJUSTED']
        if 'CDOM_ADJUSTED' in self.__floatdict__.keys():
            self.CDOM_ADJUSTED = self.__floatdict__['CDOM_ADJUSTED']

        # not naturally gridded variables
        self.CYCLE_GRID     = self.__floatdict__['CYCLE_GRID']
        self.SDN_GRID       = self.__floatdict__['SDN_GRID']
        self.LATITUDE_GRID  = self.__floatdict__['LATITUDE_GRID']
        self.LONGITUDE_GRID = self.__floatdict__['LONGITUDE_GRID']
    
    def to_dict(self):
        return self.__floatdict__.copy()
    
    def to_dataframe(self):
        import pandas as pd

        df = pd.DataFrame()
        df['CYCLE']     = self.CYCLE_GRID
        df['SDN']       = self.SDN_GRID
        df['LATITUDE']  = self.LATITUDE_GRID
        df['LONGITUDE'] = self.LONGITUDE_GRID
        df['PRES']      = self.PRES
        df['TEMP']      = self.TEMP
        df['PSAL']      = self.PSAL
        df['PDEN']      = self.PDEN
        if 'DOXY' in self.__floatdict__.keys():
            df['DOXY']      = self.DOXY
        if 'CHLA' in self.__floatdict__.keys():
            df['CHLA']      = self.CHLA
        if 'BBP700' in self.__floatdict__.keys():
            df['BBP700']    = self.BBP700
        if 'CDOM' in self.__floatdict__.keys():
            df['CDOM']      = self.CDOM
        if 'DOXY_ADJUSTED' in self.__floatdict__.keys():
            df['DOXY_ADJUSTED']      = self.DOXY_ADJUSTED
        if 'CHLA_ADJUSTED' in self.__floatdict__.keys():
            df['CHLA_ADJUSTED']      = self.CHLA_ADJUSTED
        if 'BBP700_ADJUSTED' in self.__floatdict__.keys():
            df['BBP700_ADJUSTED']    = self.BBP700_ADJUSTED
        if 'CDOM_ADJUSTED' in self.__floatdict__.keys():
            df['CDOM_ADJUSTED']      = self.CDOM_ADJUSTED

        self.df = df

        return df.copy()

    def get_track(self):
        self.track = track(self.__floatdict__)

        return self.track.copy()

    def get_ncep(self):

        if not hasattr(self, 'track'):
            self.get_track()

        self.NCEP, self.__NCEPweights__ = ncep_to_float_track('pres', self.track, local_path=self.ncep_path)
        
        return self.NCEP.copy()

    def get_woa(self):

        if not hasattr(self, 'track'):
            self.get_track()
        
        self.z_WOA, self.WOA, self.__WOAweights__ = woa_to_float_track(self.track, 'O2sat', local_path=self.woa_path)

        return self.WOA.copy()

    def calc_gains(self, ref='NCEP', zlim=25.):

        if not hasattr(self, 'track'):
            self.get_track()

        if ref == 'NCEP':
            # check if reference data is already calculated
            if not hasattr(self, 'NCEP'):
                self.get_ncep()

            ix = [c in np.unique(self.__floatdict__['TRAJ_CYCLE']) for c in self.CYCLES]
            self.NCEP_PPOX = unit.atmos_pO2(self.NCEP[ix], unit.pH2O(get_var_by('TEMP_DOXY', 'TRAJ_CYCLE', self.__floatdict__)))/100
            self.__NCEPgains__, self.__NCEPfloatref__ = calc_gain(self.__floatdict__, self.NCEP_PPOX)
            self.gains = self.__NCEPgains__

        if ref == 'WOA':
            # check if reference data is already calculated
            if not hasattr(self, 'WOA'):
                self.get_woa()

            self.__WOAgains__, self.__WOAfloatref__, self.__WOAref__ = calc_gain(self.__floatdict__, dict(z=self.z_WOA, WOA=self.WOA), inair=False, zlim=zlim)
            self.gains = self.__WOAgains__
        
        return self.gains.copy()

    def plot(self, kind, **kwargs):

        if kind == 'gain':
            ref = kwargs['ref']
            if ref == 'NCEP':
                g = fplt.gainplot(self.SDN, self.__NCEPfloatref__[:,2], self.NCEP_PPOX, self.__NCEPgains__, ref)
            elif ref == 'WOA':
                g = fplt.gainplot(self.SDN, self.__WOAfloatref__[:,2], self.__WOAref__, self.__WOAgains__, ref)

        elif kind == 'cscatter':
            var = kwargs.pop('varname')

            if not hasattr(self, 'df'):
                self.to_dataframe()

            g = fplt.var_cscatter(self.df, varname=var, **kwargs)

        elif kind == 'profiles':
            varlist = kwargs.pop('varlist')

            if not hasattr(self, 'df'):
                self.to_dataframe()

            g = fplt.profiles(self.df, varlist=varlist, **kwargs)

        return g

class profiles:

    set_dirs = set_dirs

    def __init__(self, floats, cycles=None, mission='B', mode='RD'):
        self.__argofiles__ = get_files(ARGO_PATH, floats, cycles=cycles)
        self.__floatdict__ = load_profiles(self.__argofiles__)

        # local path info
        self.argo_path = ARGO_PATH
        self.woa_path  = WOA_PATH
        self.ncep_path = NCEP_PATH

        # metadata and dimension variables
        self.floatName  = self.__floatdict__['floatName']
        self.floatType  = self.__floatdict__['floatType']
        self.N_LEVELS   = self.__floatdict__['N_LEVELS']
        self.CYCLE      = self.__floatdict__['CYCLE']
        self.CYCLE_GRID = self.__floatdict__['CYCLE_GRID']

        # time and location data
        self.SDN       = self.__floatdict__['SDN']
        self.SDN_GRID       = self.__floatdict__['SDN_GRID']
        self.LATITUDE  = self.__floatdict__['LATITUDE']
        self.LATITUDE_GRID  = self.__floatdict__['LATITUDE_GRID']
        self.LONGITUDE = self.__floatdict__['LONGITUDE']
        self.LONGITUDE_GRID = self.__floatdict__['LONGITUDE_GRID']

        self.WMO = self.__floatdict__['WMO']

        # core variables
        self.PRES    = self.__floatdict__['PRES']
        # self.PRES_QC = self.__floatdict__['PRES_QC']
        # self.TEMP    = self.__floatdict__['TEMP']
        # self.TEMP_QC = self.__floatdict__['TEMP_QC']
        # self.PSAL    = self.__floatdict__['PSAL']
        # self.PSAL_QC = self.__floatdict__['PSAL_QC']
        # potential density
        # self.PDEN = sw.pden(self.PSAL, self.TEMP, self.PRES) - 1000

        # bgc variables - not necessarily all there so check if the fields exist
        if 'DOXY' in self.__floatdict__.keys():
            self.DOXY      = self.__floatdict__['DOXY']
            self.DOXY_QC   = self.__floatdict__['DOXY_QC']
        if 'CHLA' in self.__floatdict__.keys():
            self.CHLA      = self.__floatdict__['CHLA']
            self.CHLA_QC   = self.__floatdict__['CHLA_QC']
        if 'BBP700' in self.__floatdict__.keys():
            self.BBP700    = self.__floatdict__['BBP700']
            self.BBP700_QC = self.__floatdict__['BBP700_QC']
        if 'CDOM' in self.__floatdict__.keys():
            self.CDOM      = self.__floatdict__['CDOM']
            self.CDOM_QC   = self.__floatdict__['CDOM_QC']
        
        # adjusted variables
        if 'DOXY_ADJUSTED' in self.__floatdict__.keys():
            self.DOXY_ADJUSTED      = self.__floatdict__['DOXY_ADJUSTED']
            self.DOXY_ADJUSTED_QC   = self.__floatdict__['DOXY_ADJUSTED_QC']
        if 'CHLA_ADJUSTED' in self.__floatdict__.keys():
            self.CHLA_ADJUSTED      = self.__floatdict__['CHLA_ADJUSTED']
            self.CHLA_ADJUSTED_QC   = self.__floatdict__['CHLA_ADJUSTED_QC']
        if 'BBP700_ADJUSTED' in self.__floatdict__.keys():
            self.BBP700_ADJUSTED    = self.__floatdict__['BBP700_ADJUSTED']
            self.BBP700_ADJUSTED_QC = self.__floatdict__['BBP700_ADJUSTED_QC']
        if 'CDOM_ADJUSTED' in self.__floatdict__.keys():
            self.CDOM_ADJUSTED      = self.__floatdict__['CDOM_ADJUSTED']
            self.CDOM_ADJUSTED_QC   = self.__floatdict__['CDOM_ADJUSTED_QC']
             
    def to_dict(self):
        return self.__floatdict__.copy()
    
    def to_dataframe(self):
        import pandas as pd

        df = pd.DataFrame()
        df['CYCLE']     = self.CYCLE_GRID
        df['SDN']       = self.SDN_GRID
        df['WMO']       = self.WMO
        df['LATITUDE']  = self.LATITUDE_GRID
        df['LONGITUDE'] = self.LONGITUDE_GRID
        df['PRES']      = self.PRES
        # df['TEMP']      = self.TEMP
        # df['PSAL']      = self.PSAL
        # df['PDEN']      = self.PDEN
        if 'DOXY' in self.__floatdict__.keys():
            df['DOXY']      = self.DOXY
            df['DOXY_QC']   = self.DOXY_QC
        if 'CHLA' in self.__floatdict__.keys():
            df['CHLA']      = self.CHLA
            df['CHLA_QC']   = self.CHLA_QC
        if 'BBP700' in self.__floatdict__.keys():
            df['BBP700']    = self.BBP700
            df['BBP700_QC'] = self.BBP700_QC
        if 'CDOM' in self.__floatdict__.keys():
            df['CDOM']      = self.CDOM
            df['CDOM_QC']   = self.CDOM_QC
        if 'DOXY_ADJUSTED' in self.__floatdict__.keys():
            df['DOXY_ADJUSTED']      = self.DOXY_ADJUSTED
            df['DOXY_ADJUSTED_QC']   = self.DOXY_ADJUSTED_QC
        if 'CHLA_ADJUSTED' in self.__floatdict__.keys():
            df['CHLA_ADJUSTED']      = self.CHLA_ADJUSTED
            df['CHLA_ADJUSTED_QC']   = self.CHLA_ADJUSTED_QC
        if 'BBP700_ADJUSTED' in self.__floatdict__.keys():
            df['BBP700_ADJUSTED']    = self.BBP700_ADJUSTED
            df['BBP700_ADJUSTED_QC'] = self.BBP700_ADJUSTED_QC
        if 'CDOM_ADJUSTED' in self.__floatdict__.keys():
            df['CDOM_ADJUSTED']      = self.CDOM_ADJUSTED
            df['CDOM_ADJUSTED_QC']   = self.CDOM_ADJUSTED_QC

        self.df = df

        return self.df.copy()

    def get_track(self):
        self.track = track(self.__floatdict__)
        return self.track

    def get_ncep(self):

        if not hasattr(self, 'track'):
            self.get_track()
        self.NCEP = ncep_to_float_track('pres', self.track, local_path=self.ncep_path)
        
        return self.NCEP

    def get_woa(self):

        if not hasattr(self, 'track'):
            self.get_track()
        
        self.z_WOA, self.WOA = woa_to_float_track(self.track, 'O2sat', local_path=self.woa_path)

        return self.WOA

    def calc_gains(self, ref='NCEP'):

        if not hasattr(self, 'track'):
            self.get_track()

        if ref == 'NCEP':
            sys.stdout.write('Function not built yet, returning None\n')
            self.gains = None

        if ref == 'WOA':
            # check if reference data is already calculated
            if not hasattr(self, 'WOA'):
                self.get_woa()

            self.__WOAgains__, self.__WOAfloatref__, self.__WOAref__ = calc_gain(self.__floatdict__, dict(z=self.z_WOA, WOA=self.WOA), inair=False)
            self.gains = self.__WOAgains__
        
        return self.gains

# ----------------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------------

def apply_qc_adjustment():

    return None

def get_files(local_path, wmo_numbers, cycles=None, mission='B', mode='RD'):
    local_path = Path(local_path)

    subset_index = __bgcindex__[__bgcindex__.wmo.isin(wmo_numbers)]
    if cycles is not None:
        subset_index = subset_index[subset_index.cycle.isin(cycles)]
    wcs = ['*' + a + b + '*.nc' for a in mission for b in mode]
    wcs = [w.replace('C','') for w in wcs]

    matches = [fn for sub in [fnmatch.filter(subset_index.file, w) for w in wcs] for fn in sub]
    subset_index = subset_index[subset_index.file.isin(matches)]
    local_files = [(local_path / dac / str(wmo) / 'profiles' / fn.split('/')[-1]).as_posix() for dac, wmo, fn in zip(subset_index.dac, subset_index.wmo, subset_index.file)]

    return local_files

def get_vars(files):

    nc = Dataset(Path(files[0]), 'r')
    varnames = set(nc.variables.keys())
    nc.close()

    for fn in files:
        nc = Dataset(Path(fn), 'r')
        varnames = varnames.intersection(nc.variables.keys())
        nc.close()

    varnames = list(varnames)

    return varnames

def read_qc(flags):

    decode_flags = np.array([f.decode('utf-8') for f in flags])
    decode_flags[decode_flags == ' '] = '4'

    out_flags = np.array([int(f) for f in decode_flags])

    return out_flags

def load_argo(local_path, wmo, grid=False, verbose=False):
    '''
    Function to load in all data from a single float, using BRtraj, meta,
    and Sprof files
    
    INPUT:
            local_path: local path of float data
            wmo: float ID number
    
    OUTPUT:
            floatData: python dict() object with the following fields
                floatName: WMO number, from input
                floatType: Kind of float (APEX, ARVOR, etc.)
                N_LEVELS: Number of depth levels, Argo dimension N_LEVELS
                N_CYCLES: Number of profiles, Argo dimension N_PROF
                CYCLES: Array from 1 to N_CYCLES
                LATITUDE: Latitude (-90, 90) for each profile
                LONGITUDE: Longitude (-180, 180) for each profile
                SDN: Serial Date Number for each profile
                PRES: Pressure (dbar), compressed to vector (1D array)
                TEMP: Temperature (deg C)
                PSAL: Salinity (psu)

            if the variables are available, it will also contain:
                DOXY: Dissolved Oxygen (micromole/kg)
                O2sat: Oxygen percent saturation (%)
                PPOX_DOXY: Oxygen partial pressure (mbar) [if avail.]
                TRAJ_CYCLE: Cycle number for PPOX_DOXY [if avail.]
                inair: Boolean to indicate if in-air data exists
            
            for all the variables listen above, there will also exist
                <PARAM>_QC fields for quality flags, and <PARAM>_ADJUSTED
                fields if they exist
    
                *** CYCLES, LATITUDE, LONGITUDE, AND SDN ALL ALSO HAVE ***
                ***     ANALOGOUS <VAR>_GRID FIELDS THAT MATCH THE     ***
                ***   DIMENSION OF PRES, TEMP, PSAL, DOXY, AND O2SAT   ***
    
    AUTHOR:   Christopher Gordon
              Fisheries and Oceans Canada
              chris.gordon@dfo-mpo.gc.ca
    
    ACKNOWLEDGEMENT: this code is adapted from the SOCCOM SAGE_O2Argo matlab
    code, available via https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING,
    written by Tanya Maurer & Josh Plant
    
    LAST UPDATE: 29-04-2020
    
    CHANGE LOG:
    
    22-04-2020: updated so that pressure mask determines all variables - need
    to add all quality flags to output
    
    29-04-2020: switched file/path handling from os module to pathlib
    '''

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
        if verbose:
            sys.stdout.write('Continuing without BRtraj file\n')
    elif BRtraj.exists():
        BRtraj_nc = Dataset(BRtraj, 'r')
        if 'PPOX_DOXY' not in BRtraj_nc.variables.keys():
            BRtraj_flag = False
            if verbose:
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
    floatData = dict(floatName=wmo, N_CYCLES=N, N_LEVELS=M)
    
    mask = Sprof_nc.variables['PRES'][:].mask
    mask_vars = ['TEMP','PSAL']
    if 'DOXY' in Sprof_nc.variables.keys():
        mask_vars = mask_vars + ['DOXY']
    if 'DOXY_ADJUSTED' in Sprof_nc.variables.keys():
        mask_vars = mask_vars + ['DOXY_ADJUSTED']

    for v in mask_vars:
        mask = np.logical_or(mask, Sprof_nc.variables[v][:].mask)

    if not 'CYCLE_NUMBER' in Sprof_nc.variables.keys():
        floatData['CYCLES'] = np.arange(1,N+1)
    else:
        floatData['CYCLES'] = Sprof_nc.variables['CYCLE_NUMBER'][:].data

    # load in variables that will be in every file
    floatData['PRES'] = Sprof_nc.variables['PRES'][:].data
    floatData['TEMP'] = Sprof_nc.variables['TEMP'][:].data
    floatData['PSAL'] = Sprof_nc.variables['PSAL'][:].data
    floatData['SDN']  = Sprof_nc.variables['JULD'][:].data + pl.datestr2num('1950-01-01')
    floatData['LATITUDE']  = Sprof_nc.variables['LATITUDE'][:].data
    floatData['LONGITUDE'] = Sprof_nc.variables['LONGITUDE'][:].data

    # loop through other possible BGC variables
    bgc_vars = ['DOXY', 'CHLA', 'BBP700', 'CDOM', 'NITRATE', 'DOWNWELLING_IRRADIANCE']
    for v in bgc_vars:
        if v in Sprof_nc.variables.keys():
            floatData[v] = Sprof_nc.variables[v][:].data

    for v in floatData.keys():
        v_qc = v + '_QC'
        if v_qc in Sprof_nc.variables.keys():
            floatData[v_qc] = Sprof_nc.variables[v_qc][:].data
        v_adj = v + '_ADJUSTED'
        if v_adj in Sprof_nc.variables.keys():
            floatData[v_adj] = Sprof_nc.variables[v_adj][:].data
            floatData[v_adj + '_QC'] = Sprof_nc.variabes[v_adj + '_QC'][:].data

    if grid:
        ftype = ''
        for let in meta_nc.variables['PLATFORM_TYPE'][:].compressed():
            ftype = ftype + let.decode('UTF-8')
        floatData['floatType'] = ftype

        floatData['SDN_GRID']       = np.tile(floatData['SDN'],(M,1)).T
        floatData['CYCLE_GRID']     = np.tile(floatData['CYCLES'],(M,1)).T
        floatData['LATITUDE_GRID']  = np.tile(floatData['LATITUDE'],(M,1)).T
        floatData['LONGITUDE_GRID'] = np.tile(floatData['LONGITUDE'],(M,1)).T

    floatData['O2Sat'] = 100*floatData['DOXY']/unit.oxy_sol(floatData['PSAL'], floatData['TEMP'], unit='micromole/kg')

    if BRtraj_flag:
        ppox_doxy               = BRtraj_nc.variables['PPOX_DOXY'][:]
        floatData['PPOX_DOXY']  = ppox_doxy.compressed()
        floatData['TEMP_DOXY']  = BRtraj_nc.variables['TEMP_DOXY'][:].data
        floatData['TRAJ_CYCLE'] = BRtraj_nc.variables['CYCLE_NUMBER'][:].data
        floatData['inair']      = True
    else:
        floatData['inair']      = False

    return floatData

def load_profile(fn):
    '''
    Function to load a singe Argo profile file into a dict() object
    NOTE: Deprecated, use load_profiles instead, which can handle multiple
    profile files at once, but produces the same result for just one. 

    AUTHOR:   Christopher Gordon
              Fisheries and Oceans Canada
              chris.gordon@dfo-mpo.gc.ca
    
    LAST UPDATE: 29-04-2020
    
    CHANGE LOG:
    
    22-04-2020: updated so that pressure mask determines all variables - need
    to add all quality flags to output
    
    29-04-2020: switched file/path handling from os module to pathlib

    24-06-2020: deprecated, re-wrote as load_profiles()
    '''

    # # try to load the profile as absolute path or relative path
    # try:
    #     nc = Dataset(fn, 'r')
    # except:
    #     try:
    #         nc = Dataset(Path(ARGO_PATH) / fn, 'r')
    #     except:
    #         raise FileNotFoundError('No such file {} or {}'.format(fn, str(Path(ARGO_PATH) / fn)))

    # wmo = ''
    # for let in nc.variables['PLATFORM_NUMBER'][:].compressed():
    #     wmo = wmo + let.decode('UTF-8')

    # cycle = nc.variables['CYCLE_NUMBER'][:].compressed()[0]

    # # number of profile cycles
    # M = nc.dimensions['N_LEVELS'].size
    # # beginning of output dict with basic info, following variables in SAGEO2
    # floatData = dict(floatName=wmo, N_LEVELS=M, CYCLE=cycle)

    # ftype = ''
    # for let in nc.variables['PLATFORM_TYPE'][:].compressed():
    #     ftype = ftype + let.decode('UTF-8')

    # floatData['floatType'] = ftype

    # # load in variables that will be in every file
    # floatData['PRES'] = nc.variables['PRES'][:].data
    # floatData['TEMP'] = nc.variables['TEMP'][:].data
    # floatData['PSAL'] = nc.variables['PSAL'][:].data
    # floatData['SDN']  = nc.variables['JULD'][:].data + pl.datestr2num('1950-01-01')
    # floatData['LATITUDE']  = nc.variables['LATITUDE'][:].data
    # floatData['LONGITUDE'] = nc.variables['LONGITUDE'][:].data

    # # loop through other possible BGC variables
    # bgc_vars = ['DOXY', 'CHLA', 'BBP700', 'CDOM', 'NITRATE', 'DOWNWELLING_IRRADIANCE']
    # for v in bgc_vars:
    #     if v in nc.variables.keys():
    #         floatData[v] = nc.variables[v][:].data

    # for v in floatData.keys():
    #     v_qc = v + '_QC'
    #     if v_qc in nc.variables.keys():
    #         floatData[v_qc] = nc.variables[v_qc][:].data
    #     v_adj = v + '_ADJUSTED'
    #     if v_adj in nc.variables.keys():
    #         floatData[v_adj] = nc.variables[v_adj][:].data
    #         floatData[v_adj + '_QC'] = nc.variabes[v_adj + '_QC'][:].data

    # return floatData

def load_profiles(files):

    common_variables = get_vars(files)

    floatData = dict(
        floatName=[], N_LEVELS=[], N_PROF=[], CYCLE=np.array([]), floatType=[]
    )

    for v in ['PRES', 'TEMP', 'PSAL', 'SDN']:
        floatData[v] = np.array([])
        floatData[v + '_QC'] = np.array([])
    
    floatData.pop('TEMP')
    floatData.pop('TEMP_QC')
    floatData.pop('PSAL')
    floatData.pop('PSAL_QC')
    
    for v in ['WMO', 'LATITUDE', 'LONGITUDE', 'POSITION_QC', 'SDN_GRID', 'LATITUDE_GRID', 'LONGITUDE_GRID', 'CYCLE_GRID']:
        floatData[v] = np.array([])

    for v in ['DOXY', 'CHLA', 'BBP700', 'CDOM', 'NITRATE', 'DOWNWELLING_IRRADIANCE']:
        if v in common_variables:
            floatData[v] = np.array([])
            floatData[v + '_QC'] = np.array([])
            if v + '_ADJUSTED' in common_variables:
                floatData[v + '_ADJUSTED'] = np.array([])
                floatData[v + '_ADJUSTED' + '_QC'] = np.array([])

    for fn in files:
        # try to load the profile as absolute path or relative path
        try:
            nc = Dataset(fn, 'r')
        except:
            try:
                nc = Dataset(Path(ARGO_PATH) / fn, 'r')
            except:
                raise FileNotFoundError('No such file {} or {}'.format(fn, str(Path(ARGO_PATH) / fn)))

        # number of profile cycles
        M = nc.dimensions['N_LEVELS'].size
        N = nc.dimensions['N_PROF'].size

        wmo = ''
        if N > 1:
            for let in nc.variables['PLATFORM_NUMBER'][:][0,:].compressed():
                wmo = wmo + let.decode('UTF-8')
        else:
            for let in nc.variables['PLATFORM_NUMBER'][:].compressed():
                wmo = wmo + let.decode('UTF-8')

        cycle = nc.variables['CYCLE_NUMBER'][:].data.flatten()

        ftype = ''
        for let in nc.variables['PLATFORM_TYPE'][:].compressed():
            ftype = ftype + let.decode('UTF-8')

        floatData['floatName']  = floatData['floatName'] + [int(wmo)]
        floatData['N_LEVELS']   = floatData['N_LEVELS']  + [M]
        floatData['N_PROF']     = floatData['N_PROF']    + [N]
        floatData['CYCLE']      = np.append(floatData['CYCLE'], cycle)
        floatData['CYCLE_GRID'] = np.append(floatData['CYCLE_GRID'], np.array(N*M*[cycle[0]]))
        floatData['floatType']  = floatData['floatType'] + [ftype]
        floatData['WMO']        = np.append(floatData['WMO'], np.array(M*N*[wmo]))

        # load in variables that will be in every file
        floatData['PRES'] = np.append(floatData['PRES'], nc.variables['PRES'][:].data.flatten())
        # floatData['TEMP'] = np.append(floatData['TEMP'], nc.variables['TEMP'][:].data.flatten())
        # floatData['PSAL'] = np.append(floatData['PSAL'], nc.variables['PSAL'][:].data.flatten())
        floatData['SDN']  = np.append(floatData['SDN'], nc.variables['JULD'][:].data.flatten() + pl.datestr2num('1950-01-01'))
        floatData['SDN_QC'] = np.append(floatData['SDN_QC'], read_qc(nc.variables['JULD_QC'][:].data.flatten()))
        floatData['SDN_GRID'] = np.append(floatData['SDN_GRID'], np.array(N*M*[np.nanmean(nc.variables['JULD'][:].data.flatten() + pl.datestr2num('1950-01-01'))]))
        floatData['LATITUDE'] = np.append(floatData['LATITUDE'], nc.variables['LATITUDE'][:].data.flatten())
        floatData['LATITUDE_GRID'] = np.append(floatData['LATITUDE_GRID'], np.array(N*M*[np.nanmean(nc.variables['LATITUDE'][:].data.flatten())]))
        floatData['LONGITUDE'] = np.append(floatData['LONGITUDE'], nc.variables['LONGITUDE'][:].data.flatten())
        floatData['LONGITUDE_GRID'] = np.append(floatData['LONGITUDE_GRID'], np.array(N*M*[np.nanmean(nc.variables['LONGITUDE'][:].data.flatten())]))
        floatData['POSITION_QC'] = np.append(floatData['POSITION_QC'], read_qc(nc.variables['POSITION_QC'][:].data.flatten()))

        # loop through other possible BGC variables
        bgc_vars = ['DOXY', 'CHLA', 'BBP700', 'CDOM', 'NITRATE', 'DOWNWELLING_IRRADIANCE']
        for v in bgc_vars:
            if v in common_variables:
                floatData[v] = np.append(floatData[v], nc.variables[v][:].data)
            v_adj = v + '_ADJUSTED'
            if v_adj in common_variables:
                floatData[v_adj] = np.append(floatData[v_adj], nc.variables[v_adj][:].data.flatten())

        for v in floatData.keys():
            v_qc = v + '_QC'
            if v_qc in common_variables:
                floatData[v_qc] = np.append(floatData[v_qc], read_qc(nc.variables[v_qc][:].data.flatten()))

    return floatData

def clean(float_data):

    clean_float_data = float_data.copy()

    for qc_key in any('_QC' for key in clean_float_data.keys()):
        data_key   = qc_key.replace('_QC','')
        good_index = np.logical_or(np.logical_or(clean_float_data[qc_key] < 3 | clean_float_data[qc_key] == 5), clean_float_data[qc_key] == 8)
        bad_index  = np.invert(good_index)

        clean_float_data[data_key][bad_index] = np.nan

    return clean_float_data

def track(float_data):
    # make 'track' array with columns (time, lat, lon) to be used in interpolation
    track = np.array([float_data['SDN'], float_data['LATITUDE'], float_data['LONGITUDE']]).T

    return track

def woa_to_float_track(track, param, zlim=(0,1000), local_path='./'):
    '''
    Function to load WOA18 climatological data for comparison with autonomous
    floats. Data to be interpolated along the provided track (t, lat, lon).
    Combines function load_woa_data() and interp_woa_data() for convenience,
    see documentation for those funcions for more detail.
    
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
              z: WOA depth array
              woa_interp: 2D array of requested WOA parameter (depth x time)
    
    AUTHOR:   Christopher Gordon
              Fisheries and Oceans Canada
              chris.gordon@dfo-mpo.gc.ca
    
    LAST UPDATE: 23-04-2020
    
    CHANGE LOG:
    '''

    xtrack, woa_track, woa_data = io.load_woa_data(track, param, zlim=zlim, local_path=local_path)
    woa_interp, wt, yrday = interp.interp_woa_data(xtrack, woa_track, woa_data)
    z = woa_track[0]

    return z, woa_interp, wt

def ncep_to_float_track(varname, track, local_path='./'):
    '''
    Function to load NCEP reanalysis data for comparison with autonomous
    floats. Data to be interpolated along the provided track (t, lat, lon).
    Combines function load_ncep_data() and interp_ncep_data() for convenience,
    see documentation for those funcions for more detail.
    
    INPUT:
              varname: either 'pres' (pressure) or 'rhum' (relative humidity)
              track: array with the columns (SDN, lat, lon)
    
    OUTPUT:
              z: WOA depth array
              woa_interp: 2D array of requested WOA parameter (depth x time)
    
    AUTHOR:   Christopher Gordon
              Fisheries and Oceans Canada
              chris.gordon@dfo-mpo.gc.ca
    
    LAST UPDATE: 29-04-2020
    
    CHANGE LOG:
    '''

    xtrack, ncep_track, data = io.load_ncep_data(track, varname, local_path=local_path)
    ncep_interp, wt = interp.interp_ncep_data(xtrack, ncep_track, data)

    return ncep_interp, wt


def calc_gain(data, ref, inair=True, zlim=25., verbose=False):
    '''
    Calculate the gain for each profile by comparing float oxygen data to a
    reference data set, either NCEP for in-air or WOA surface data if in-air
    comparison is not available.
    
    INPUT:
              data: float data dict object, output from load_argo_data()
              ref: reference data set, either NCEP pO2 or WOA O2sat
              inair: boolean flag to indicate if comparison to NCEP in-air
                  data or WOA surface data should be done, default to
                  in-air, but function also performs check
              zlim: lower limit to define as 'surface' and take mean within,
                    default value 25 dbar, for use only when inair is False
    
    OUTPUT:
              g: vector of gains
              surf_data: array of float surface stats (cycle, N, mean, std)
    
    AUTHOR:   Christopher Gordon
              Fisheries and Oceans Canada
              chris.gordon@dfo-mpo.gc.ca
    
    LAST UPDATE: 23-04-2020
    
    CHANGE LOG:
    '''

    # check which reference data to use
    if inair and 'PPOX_DOXY' not in data.keys():
        raise ValueError('Flag ''inair'' set to True but partial pressure data not available')

    if inair:
        if verbose:
            sys.stdout.write('\nCalculating gains using NCEP surface pressure and float in-air measurements...\n')
        g = np.nan*np.ones((ref.shape[0],))

        # float partial pressure measurements at each cycle
        ppox  = data['PPOX_DOXY']
        cycle = data['CYCLES']
        inair_cycle = data['TRAJ_CYCLE']

        mean_float_data = np.nan*np.ones((ref.shape[0],4))
        for i,c in enumerate(np.unique(inair_cycle)):
            subset_ppox = ppox[inair_cycle == c]
            mean_float_data[i,0] = c
            mean_float_data[i,1] = np.sum(~np.isnan(subset_ppox))
            mean_float_data[i,2] = np.nanmean(subset_ppox)
            mean_float_data[i,3] = np.nanstd(subset_ppox)

            g[i] = ref[i]/mean_float_data[i,2]

        g[g == 0] = np.nan

        return g, mean_float_data

    else:
        if verbose:
            sys.stdout.write('\nCalculating gains using WOA surface data and float O2 percent saturation...\n')
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
        
        g[g == 0] = np.nan

        return g, mean_float_data, woa_surf

def aic(data,resid):
    '''
    Function to calculate the Akiake Information Criteria (AIC) as a metric
    for assessing the appropriate number of breakpoints in the calculation of
    drifts in O2 gains.
    
    INPUT:
    
    OUTPUT:
    
    AUTHOR:   Christopher Gordon
              Fisheries and Oceans Canada
              chris.gordon@dfo-mpo.gc.ca
    
    ACKNOWLEDGEMENT: this code is adapted from the SOCCOM SAGE_O2Argo matlab
    code, available via https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING,
    written by Tanya Maurer & Josh Plant
    
    LAST UPDATE: 20-04-2020
    
    CHANGE LOG:
    '''

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
    '''
    Function to calculate the Bayesian Information Criteria (BIC) as a metric
    for assessing the appropriate number of breakpoints in the calculation of
    drifts in O2 gains.
    
    INPUT:
    
    OUTPUT:
    
    AUTHOR:   Christopher Gordon
              Fisheries and Oceans Canada
              chris.gordon@dfo-mpo.gc.ca
    
    ACKNOWLEDGEMENT: this code is adapted from the SOCCOM SAGE_O2Argo matlab
    code, available via https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING,
    written by Tanya Maurer & Josh Plant
    
    LAST UPDATE: 20-04-2020
    
    CHANGE LOG:
    '''

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
    '''
    Function to calculate the mean of one variable (v1) indexed by a second
    variable (v2) in a float_data dictionary (output of sagepy.argo), though
    it would work with any python dict
    
    INPUT:
              v1: string input of a key in float_data
              v2: string input of a key in float_data
              float_data: python dict() object
    
    OUTPUT:
              out_array: 1D numpy array with mean values
    
    AUTHOR:   Christopher Gordon
              Fisheries and Oceans Canada
              chris.gordon@dfo-mpo.gc.ca
    
    ACKNOWLEDGEMENT: this code is adapted from the SOCCOM SAGE_O2Argo matlab
    code, available via https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING,
    written by Tanya Maurer & Josh Plant
    
    LAST UPDATE: 20-04-2020
    
    CHANGE LOG:
    '''

    
    index = np.unique(float_data[v2])
    out_array = np.nan*np.ones((len(index)))
    for i,v in enumerate(index):
        out_array[i] = np.nanmean(float_data[v1][float_data[v2] == v])
    
    return out_array

def oxy_b(dt, tau):
    inv_b = 1 + 2*(tau/dt)
    return 1/inv_b

def oxy_a(dt, tau):
    return 1 - 2*oxy_b(dt, tau)

def correct_response_time(t, DO, T, thickness):

    # convert time to seconds
    t_sec = t*24*60*60

    # array for the loop
    N = DO.shape[0]
    mean_oxy  = np.array((N-1)*[np.nan])
    mean_time = t_sec[:-1] + np.diff(t_sec)/2
    mean_temp = T[:-1] + np.diff(T)/2

    # load temperature, boundary layer thickness, and tau matrix from 
    # look-up table provided in the supplement to Bittig and Kortzinger (2017)
    lut_data = np.loadtxt(Path('LUT/T_lL_tau_3830_4330.dat'))
    lut_lL = lut_data[0,1:]
    lut_T  = lut_data[1:,0]
    tau100 = lut_data[1:,1:]
    thickness = thickness*np.ones((N-1,))

    # translate boundary layer thickness to temperature dependent tau
    f_thickness = RectBivariateSpline(lut_T, lut_lL, tau100, kx=1, ky=1)
    tau_T = np.squeeze(f_thickness(thickness[0], mean_temp))

    # loop through oxygen data
    for i in range(N-1):
        dt = t_sec[i+1] - t_sec[i]

        # do the correction using the mean filter, get the mean time
        mean_oxy[i]  = (1/(2*oxy_b(dt, tau_T[i])))*(DO[i+1] - oxy_a(dt, tau_T[i])*DO[i])
    
    # interpolate back to original times for output
    f = interp1d(mean_time, mean_oxy, kind='linear', bounds_error=False, fill_falue='extrapolate')
    DO_out = f(t_sec)

    return DO_out

def correct_response_time_Tconst(t, DO, tau):

    # array for the loop
    N = DO.shape[0]
    mean_oxy  = np.array((N-1)*[np.nan])
    mean_time = np.array((N-1)*[np.nan])

    # convert time to seconds
    t_sec = t*24*60*60

    # loop through oxygen data
    for i in range(N-1):
        dt = t_sec[i+1] - t_sec[i]

        # do the correction using the mean filter, get the mean time
        mean_oxy[i]  = (1/(2*oxy_b(dt, tau)))*(DO[i+1] - oxy_a(dt, tau)*DO[i])
        mean_time[i] = t_sec[i] + dt/2
    
    # interpolate back to original times for output
    f = interp1d(mean_time, mean_oxy, kind='linear', bounds_error=False, fill_falue='extrapolate')
    DO_out = f(t_sec)

    return DO_out