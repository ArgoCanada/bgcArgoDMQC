#!/usr/bin/python

import sys
import copy
import warnings
from pathlib import Path
import fnmatch
import time

import numpy as np
from scipy.interpolate import interp1d, RectBivariateSpline

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# soft attempt to load gsw, but allow for seawater as well
try: 
    import gsw
    flagSA = True
except:
    try:
        # if this also fails, just load gsw to throw the error
        from seawater import pden
        flagSA = False
        warnings.warn('gsw package for thermodynamic equations of seawater not installed, attempting to load seawater package, however seawater is deprecated in favour of gsw-python, see https://teos-10.github.io/GSW-Python/\n')
    except:
        import gsw

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
global REF_PATH
REF_PATH = Path(__file__).parent.absolute() / 'ref'

def set_dirs(argo_path=ARGO_PATH, woa_path=WOA_PATH, ncep_path=NCEP_PATH):

    global ARGO_PATH
    ARGO_PATH = argo_path
    global WOA_PATH
    WOA_PATH  = woa_path
    global NCEP_PATH
    NCEP_PATH = ncep_path

def get_index(index='bgc'):
    if index == 'bgc':
        return __bgcindex__
    elif index == 'global':
        if '__globalindex__' not in globals():
            global __globalindex__
            __globalindex__ = io.read_index(mission='C')
        return __globalindex__
    elif index == 'synthetic':
        if '__synthindex__' not in globals():
            global __synthindex__
            __synthindex__ = io.read_index(mission='S')
        return __synthindex__

def get_dac(wmo):

    if '__globalindex__' not in globals():
            global __globalindex__
            __globalindex__ = io.read_index(mission='C')
    
    dac = __globalindex__[__globalindex__.wmo == wmo].dac.iloc[0]

    return dac

# ----------------------------------------------------------------------------
# FLOAT CLASS
# ----------------------------------------------------------------------------

class sprof:

    set_dirs = set_dirs

    def __init__(self, wmo, keep_fillvalue=False):
        self.__floatdict__, self.__Sprof__, self.__BRtraj__, self.__meta__ = load_argo(ARGO_PATH, wmo, grid=True)
        self.__rawfloatdict__ = self.__floatdict__

        # local path info
        self.argo_path = ARGO_PATH
        self.woa_path  = WOA_PATH
        self.ncep_path = NCEP_PATH

        self.assign(self.__floatdict__)
        if not keep_fillvalue:
            self.rm_fillvalue()

    def assign(self, floatdict):

        # metadata and dimension variables
        self.floatName  = floatdict['floatName']
        self.floatType  = floatdict['floatType']
        self.N_CYCLES   = floatdict['N_CYCLES']
        self.N_LEVELS   = floatdict['N_LEVELS']
        self.CYCLE      = floatdict['CYCLES']
        self.CYCLE_GRID = floatdict['CYCLE_GRID']

        # time and location data
        self.SDN       = floatdict['SDN']
        self.SDN_GRID  = floatdict['SDN_GRID']
        self.LATITUDE  = floatdict['LATITUDE']
        self.LATITUDE_GRID = floatdict['LATITUDE_GRID']
        self.LONGITUDE = floatdict['LONGITUDE']
        self.LONGITUDE_GRID = floatdict['LONGITUDE_GRID']

        self.WMO = floatdict['WMO']

        # core variables
        self.PRES    = floatdict['PRES']
        self.PRES_QC = floatdict['PRES_QC']
        self.TEMP    = floatdict['TEMP']
        self.TEMP_QC = floatdict['TEMP_QC']
        self.PSAL    = floatdict['PSAL']
        self.PSAL_QC = floatdict['PSAL_QC']
        # potential density
        if flagSA:
            self.PDEN = gsw.pot_rho_t_exact(gsw.SA_from_SP(self.PSAL, self.PRES, self.LONGITUDE_GRID, self.LATITUDE_GRID), self.TEMP, self.LONGITUDE_GRID, self.LATITUDE_GRID) - 1000
        else:
            self.PDEN = pden(self.PSAL, self.TEMP, self.PRES, 0) - 1000

        # bgc variables - not necessarily all there so check if the fields exist
        if 'DOXY' in floatdict.keys():
            self.DOXY      = floatdict['DOXY']
            self.DOXY_QC   = floatdict['DOXY_QC']
        if 'CHLA' in floatdict.keys():
            self.CHLA      = floatdict['CHLA']
            self.CHLA_QC   = floatdict['CHLA_QC']
        if 'BBP700' in floatdict.keys():
            self.BBP700    = floatdict['BBP700']
            self.BBP700_QC = floatdict['BBP700_QC']
        if 'CDOM' in floatdict.keys():
            self.CDOM      = floatdict['CDOM']
            self.CDOM_QC   = floatdict['CDOM_QC']
        
        # adjusted variables
        if 'DOXY_ADJUSTED' in floatdict.keys():
            self.DOXY_ADJUSTED      = floatdict['DOXY_ADJUSTED']
            self.DOXY_ADJUSTED_QC   = floatdict['DOXY_ADJUSTED_QC']
        if 'CHLA_ADJUSTED' in floatdict.keys():
            self.CHLA_ADJUSTED      = floatdict['CHLA_ADJUSTED']
            self.CHLA_ADJUSTED_QC   = floatdict['CHLA_ADJUSTED_QC']
        if 'BBP700_ADJUSTED' in floatdict.keys():
            self.BBP700_ADJUSTED    = floatdict['BBP700_ADJUSTED']
            self.BBP700_ADJUSTED_QC = floatdict['BBP700_ADJUSTED_QC']
        if 'CDOM_ADJUSTED' in floatdict.keys():
            self.CDOM_ADJUSTED      = floatdict['CDOM_ADJUSTED']
            self.CDOM_ADJUSTED_QC   = floatdict['CDOM_ADJUSTED_QC']

        if 'O2Sat' in floatdict.keys():
            self.O2Sat = floatdict['O2Sat']
            self.O2Sat_QC = floatdict['O2Sat_QC']

    def rm_fillvalue(self):
        self.__nofillvaluefloatdict__ = dict_fillvalue_clean(self.__rawfloatdict__)
        self.__floatdict__ = copy.deepcopy(self.__nofillvaluefloatdict__)
        self.assign(self.__nofillvaluefloatdict__)

    def clean(self, bad_flags=None):
        self.__cleanfloatdict__ = dict_clean(self.__floatdict__, bad_flags=bad_flags)
        self.__floatdict__ = copy.deepcopy(self.__cleanfloatdict__)
        self.assign(self.__cleanfloatdict__)

    def reset(self):
        self.__floatdict__ = copy.deepcopy(self.__rawfloatdict__)
        self.assign(self.__rawfloatdict__)

    def check_doxy_range(self):
        self.__doxyrangedict__ = doxy_range_check(self.__floatdict__)
        self.__floatdict__ = self.__doxyrangedict__
        self.assign(self.__doxyrangedict__)
    
    def to_dict(self):
        return copy.deepcopy(self.__floatdict__)
    
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
        if 'O2Sat' in self.__floatdict__.keys():
            df['O2Sat']      = self.O2Sat

        self.df = df

        return copy.deepcopy(self.df)

    def get_track(self):
        self.track = track(self.__floatdict__)

        return copy.deepcopy(self.track)

    def get_ncep(self):

        if not hasattr(self, 'track'):
            self.get_track()

        self.NCEP, self.__NCEPweights__ = ncep_to_float_track('pres', self.track, local_path=self.ncep_path)
        
        return copy.deepcopy(self.NCEP)

    def get_woa(self):

        if not hasattr(self, 'track'):
            self.get_track()
        
        self.z_WOA, self.WOA, self.__WOAweights__ = woa_to_float_track(self.track, 'O2sat', local_path=self.woa_path)

        return copy.deepcopy(self.WOA)

    def calc_gains(self, ref='NCEP', zlim=25.):

        if not hasattr(self, 'track'):
            self.get_track()

        if ref == 'NCEP':
            # check if reference data is already calculated
            if not hasattr(self, 'NCEP'):
                self.get_ncep()

            pH2O = unit.pH2O(get_var_by('TEMP_DOXY', 'TRAJ_CYCLE', self.__floatdict__))

            common_cycles, c1, c2 = np.intersect1d(self.CYCLE, np.unique(self.__floatdict__['TRAJ_CYCLE']), assume_unique=True, return_indices=True)

            self.NCEP_PPOX = unit.atmos_pO2(self.NCEP[c1], pH2O[c2])/100
            self.__NCEPgains__, self.__NCEPfloatref__ = calc_gain(self.__floatdict__, self.NCEP_PPOX)
            self.gains = self.__NCEPgains__

        if ref == 'WOA':
            # check if reference data is already calculated
            if not hasattr(self, 'WOA'):
                self.get_woa()

            self.__WOAgains__, self.__WOAfloatref__, self.__WOAref__ = calc_gain(self.__floatdict__, dict(z=self.z_WOA, WOA=self.WOA), inair=False, zlim=zlim)
            self.gains = self.__WOAgains__
        
        return copy.deepcopy(self.gains)

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


    def describe(self):

        if not hasattr(self, 'df'):
            self.to_dataframe()    
        
        print(self.df.describe())

class profiles:

    set_dirs = set_dirs

    def __init__(self, floats, cycles=None, mission='B', mode='RD', keep_fillvalue=False):
        if type(floats) is int:
            floats = [floats]

        self.__argofiles__ = get_files(ARGO_PATH, floats, cycles=cycles, mission=mission, mode=mode)
        self.__floatdict__ = load_profiles(self.__argofiles__)
        self.__rawfloatdict__ = self.__floatdict__

        # local path info
        self.argo_path = ARGO_PATH
        self.woa_path  = WOA_PATH
        self.ncep_path = NCEP_PATH

        self.assign(self.__floatdict__)
        if not keep_fillvalue:
            self.rm_fillvalue()

    def assign(self, floatdict):

        # metadata and dimension variables
        self.floatName  = floatdict['floatName']
        self.floatType  = floatdict['floatType']
        self.N_LEVELS   = floatdict['N_LEVELS']
        self.CYCLE      = floatdict['CYCLES']
        self.CYCLE_GRID = floatdict['CYCLE_GRID']

        # time and location data
        self.SDN       = floatdict['SDN']
        self.SDN_GRID  = floatdict['SDN_GRID']
        self.LATITUDE  = floatdict['LATITUDE']
        self.LATITUDE_GRID  = floatdict['LATITUDE_GRID']
        self.LONGITUDE = floatdict['LONGITUDE']
        self.LONGITUDE_GRID = floatdict['LONGITUDE_GRID']

        self.WMO = floatdict['WMO']

        # core variables
        self.PRES    = floatdict['PRES']
        # self.PRES_QC = floatdict['PRES_QC']
        if 'TEMP' in floatdict.keys():
            self.TEMP    = floatdict['TEMP']
            self.TEMP_QC = floatdict['TEMP_QC']
            self.PSAL    = floatdict['PSAL']
            self.PSAL_QC = floatdict['PSAL_QC']
            # potential density
            if flagSA:
                self.PDEN = gsw.pot_rho_t_exact(gsw.SA_from_SP(self.PSAL, self.PRES, self.LONGITUDE_GRID, self.LATITUDE_GRID), self.TEMP, self.LONGITUDE_GRID, self.LATITUDE_GRID) - 1000
            else:
                self.PDEN = pden(self.PSAL, self.TEMP, self.PRES, 0) - 1000

        # bgc variables - not necessarily all there so check if the fields exist
        if 'DOXY' in floatdict.keys():
            self.DOXY      = floatdict['DOXY']
            self.DOXY_QC   = floatdict['DOXY_QC']
        if 'CHLA' in floatdict.keys():
            self.CHLA      = floatdict['CHLA']
            self.CHLA_QC   = floatdict['CHLA_QC']
        if 'BBP700' in floatdict.keys():
            self.BBP700    = floatdict['BBP700']
            self.BBP700_QC = floatdict['BBP700_QC']
        if 'CDOM' in floatdict.keys():
            self.CDOM      = floatdict['CDOM']
            self.CDOM_QC   = floatdict['CDOM_QC']
        
        # adjusted variables
        if 'DOXY_ADJUSTED' in floatdict.keys():
            self.DOXY_ADJUSTED      = floatdict['DOXY_ADJUSTED']
            self.DOXY_ADJUSTED_QC   = floatdict['DOXY_ADJUSTED_QC']
        if 'CHLA_ADJUSTED' in floatdict.keys():
            self.CHLA_ADJUSTED      = floatdict['CHLA_ADJUSTED']
            self.CHLA_ADJUSTED_QC   = floatdict['CHLA_ADJUSTED_QC']
        if 'BBP700_ADJUSTED' in floatdict.keys():
            self.BBP700_ADJUSTED    = floatdict['BBP700_ADJUSTED']
            self.BBP700_ADJUSTED_QC = floatdict['BBP700_ADJUSTED_QC']
        if 'CDOM_ADJUSTED' in floatdict.keys():
            self.CDOM_ADJUSTED      = floatdict['CDOM_ADJUSTED']
            self.CDOM_ADJUSTED_QC   = floatdict['CDOM_ADJUSTED_QC']

        if 'O2Sat' in floatdict.keys():
            self.O2Sat = floatdict['O2Sat']
            self.O2Sat_QC = floatdict['O2Sat_QC']

    def rm_fillvalue(self):
        self.__nofillvaluefloatdict__ = dict_fillvalue_clean(self.__rawfloatdict__)
        self.__floatdict__ = self.__nofillvaluefloatdict__
        self.assign(self.__nofillvaluefloatdict__)

    def clean(self, bad_flags=None):
        self.__cleanfloatdict__ = dict_clean(self.__floatdict__, bad_flags=bad_flags)
        self.__floatdict__ = self.__cleanfloatdict__
        self.assign(self.__cleanfloatdict__)

    def reset(self):
        self.__floatdict__ = self.__rawfloatdict__
        self.assign(self.__rawfloatdict__)

    def check_doxy_range(self):
        self.__doxyrangedict__ = doxy_range_check(self.__floatdict__)
        self.__floatdict__ = self.__doxyrangedict__
        self.assign(self.__doxyrangedict__)
             
    def to_dict(self):
        return copy.deepcopy(self.__floatdict__)
    
    def to_dataframe(self):
        import pandas as pd

        df = pd.DataFrame()
        df['CYCLE']     = self.CYCLE_GRID
        df['SDN']       = self.SDN_GRID
        df['WMO']       = self.WMO
        df['LATITUDE']  = self.LATITUDE_GRID
        df['LONGITUDE'] = self.LONGITUDE_GRID
        df['PRES']      = self.PRES
        df['TEMP']      = self.TEMP
        df['TEMP_QC']   = self.TEMP_QC
        df['PSAL']      = self.PSAL
        df['PSAL_QC']   = self.PSAL_QC
        df['PDEN']      = self.PDEN
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
        if 'O2Sat' in self.__floatdict__.keys():
            df['O2Sat']      = self.O2Sat
            df['O2Sat_QC']   = self.O2Sat_QC

        self.df = df

        return copy.deepcopy(self.df)

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
        
        self.z_WOA, self.WOA, self.__WOAweights__  = woa_to_float_track(self.track, 'O2sat', local_path=self.woa_path)

        return self.WOA

    def calc_gains(self, ref='WOA'):

        if not hasattr(self, 'track'):
            self.get_track()

        if ref == 'NCEP':
            sys.stdout.write('In-air data contained in BRtraj file, NCEP not a valid reference for individual profile files, returning None\n')
            self.gains = None

        if ref == 'WOA':
            # check if reference data is already calculated
            if not hasattr(self, 'WOA'):
                self.get_woa()

            self.__WOAgains__, self.__WOAfloatref__, self.__WOAref__ = calc_gain(self.__floatdict__, dict(z=self.z_WOA, WOA=self.WOA), inair=False)
            self.gains = self.__WOAgains__
        
        return self.gains

    def adjust_oxygen(self, G, eG):

        self.I_DOXY_ADJUSTED = apply_gain(self.DOXY, G)
        self.I_DOXY_ADJUSTED_QC = self.DOXY_QC
        self.I_DOXY_ADJUSTED_ERROR = calc_doxy_error(self.DOXY, G, eG)

        # metadata that needs to be recorded but I don't know where to put it yet
        SCIENTIFIC_CALIB_DATE = time.strftime('%d %b %Y')
        SCIENTIFIC_CALIB_COMMENT = 'No additional comment'
        SCIENTIFIC_CALIB_EQUATION = 'No additional equation'
        SCIENTIFIC_CALIB_COEFFICIENT = 'No additional coefficient'

        return self.I_DOXY_ADJUSTED

    def reassign_flags(self):

        return

    def assess_profile_flags(self):

        return

    def describe(self):

        if not hasattr(self, 'df'):
            self.to_dataframe()
        
        print(self.df.describe())
# ----------------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------------

def apply_gain(DOXY, G):

    DOXY_ADJUSTED = G*DOXY

    return DOXY_ADJUSTED

def calc_doxy_error(DOXY, G, eG):

    return 1

def get_files(local_path, wmo_numbers, cycles=None, mission='B', mode='RD', verbose=True):
    local_path = Path(local_path)

    if mission == 'B':
        subset_index = __bgcindex__[__bgcindex__.wmo.isin(wmo_numbers)]
    if mission == 'C':
        __coreindex__ = io.read_index(mission='C')
        subset_index = __coreindex__[__coreindex__.wmo.isin(wmo_numbers)]
    if cycles is not None:
        subset_index = subset_index[subset_index.cycle.isin(cycles)]
    wcs = ['*' + a + b + '*.nc' for a in mission for b in mode]
    wcs = [w.replace('C','') for w in wcs]

    matches = [fn for sub in [fnmatch.filter(subset_index.file, w) for w in wcs] for fn in sub]
    subset_index = subset_index[subset_index.file.isin(matches)]
    local_files = [(local_path / dac / str(wmo) / 'profiles' / fn.split('/')[-1]).as_posix() for dac, wmo, fn in zip(subset_index.dac, subset_index.wmo, subset_index.file)]

    remove_ix = []
    for i,fn in enumerate(local_files):
        if not Path(fn).exists():
            if verbose:
                sys.stdout.write('File {} does not exists locally - removing from returned list, suggest the user downloads using bgcArgo.io.get_argo(...)\n'.format(fn))
            remove_ix.append(i)
    
    if len(remove_ix) > 0:
        for ix in remove_ix[::-1]:
            local_files.pop(ix)

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

def get_worst_flag(*args):
    out_flags = np.ones(args[0].shape)

    ### block of code to find the worst flags
    if len(args) == 1:
        out_flags = args[0]
    else:
        # make an array where all data marked as good
        out_flags = np.ones(args[0].shape)
        # loop through input flags
        for flags in args:
            # loop through each datapoint flag
            for i,f in enumerate(flags):
                if f > out_flags[i] and f <= 4:
                    out_flags[i] = f
                if f in [5,8] and out_flags[i] == 1:
                    out_flags[i] = f
                if f == 4 and out_flags[i] in [5,8]:
                    out_flags[i] = f
                if f == 9:
                    out_flags[i] = 9
        
    return out_flags

def load_argo(local_path, wmo, grid=False, verbose=True):
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
    dac = get_dac(wmo)

    if type(wmo) is not str:
        wmo = str(wmo)

    # check that necessary files exist - can continue without BRtraj file but
    # need Sprof and meta files
    BRtraj = local_path / dac / wmo / '{}_BRtraj.nc'.format(wmo)
    Sprof  = local_path / dac / wmo / '{}_Sprof.nc'.format(wmo)
    meta   = local_path / dac /wmo / '{}_meta.nc'.format(wmo)

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
    floatData = dict(floatName=wmo, N_CYCLES=N, N_LEVELS=M, WMO=int(wmo))
    
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
        floatData['CYCLES'] = Sprof_nc.variables['CYCLE_NUMBER'][:].data.flatten()

    # load in variables that will be in every file
    floatData['PRES'] = Sprof_nc.variables['PRES'][:].data.flatten()
    floatData['TEMP'] = Sprof_nc.variables['TEMP'][:].data.flatten()
    floatData['PSAL'] = Sprof_nc.variables['PSAL'][:].data.flatten()
    floatData['SDN']  = Sprof_nc.variables['JULD'][:].data.flatten() + mdates.datestr2num('1950-01-01')
    floatData['LATITUDE']  = Sprof_nc.variables['LATITUDE'][:].data.flatten()
    floatData['LONGITUDE'] = Sprof_nc.variables['LONGITUDE'][:].data.flatten()

    # loop through other possible BGC variables
    bgc_vars = ['DOXY', 'CHLA', 'BBP700', 'CDOM', 'NITRATE', 'DOWNWELLING_IRRADIANCE']
    core_vars = ['PRES', 'TEMP', 'PSAL', 'POSITION']
    for v in bgc_vars:
        if v in Sprof_nc.variables.keys():
            floatData[v] = Sprof_nc.variables[v][:].data.flatten()

    for v in bgc_vars + core_vars:
        v_qc = v + '_QC'
        if v_qc in Sprof_nc.variables.keys():
            floatData[v_qc] = read_qc(Sprof_nc.variables[v_qc][:].data.flatten())
        v_adj = v + '_ADJUSTED'
        if v_adj in Sprof_nc.variables.keys():
            floatData[v_adj] = Sprof_nc.variables[v_adj][:].data.flatten()
            v_adj_qc = v_adj + '_QC'
            if v_adj_qc in Sprof_nc.variables.keys():
                floatData[v_adj_qc] = read_qc(Sprof_nc.variables[v_adj_qc][:].data.flatten())

    if grid:
        ftype = ''
        if 'PLATFORM_TYPE' in meta_nc.variables.keys():
            for let in meta_nc.variables['PLATFORM_TYPE'][:].compressed():
                ftype = ftype + let.decode('UTF-8')
        floatData['floatType'] = ftype

        floatData['SDN_GRID']       = np.tile(floatData['SDN'],(M,1)).T.flatten()
        floatData['CYCLE_GRID']     = np.tile(floatData['CYCLES'],(M,1)).T.flatten()
        floatData['LATITUDE_GRID']  = np.tile(floatData['LATITUDE'],(M,1)).T.flatten()
        floatData['LONGITUDE_GRID'] = np.tile(floatData['LONGITUDE'],(M,1)).T.flatten()

    floatData['O2Sat'] = 100*floatData['DOXY']/unit.oxy_sol(floatData['PSAL'], floatData['TEMP'], unit='micromole/kg')
    # match the fill values
    ix = np.logical_or(np.logical_or(floatData['PSAL'] >= 99999., floatData['TEMP'] >= 99999.), floatData['DOXY'] >= 99999.)
    floatData['O2Sat'][ix] = 99999.
    # get the worst QC flag from each quantity that goes into the calculation
    floatData['O2Sat_QC'] = get_worst_flag(floatData['TEMP_QC'], floatData['PSAL_QC'], floatData['DOXY_QC'])

    if BRtraj_flag:
        floatData['PPOX_DOXY']  = BRtraj_nc.variables['PPOX_DOXY'][:].data.flatten()
        floatData['TEMP_DOXY']  = BRtraj_nc.variables['TEMP_DOXY'][:].data.flatten()
        floatData['TRAJ_CYCLE'] = BRtraj_nc.variables['CYCLE_NUMBER'][:].data.flatten()
        floatData['inair']      = True
    else:
        floatData['inair']      = False

    return floatData, Sprof, BRtraj, meta

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
    # floatData['SDN']  = nc.variables['JULD'][:].data + mdates.datestr2num('1950-01-01')
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
    core_files = [fn.replace('B','') for fn in files]

    floatData = dict(
        floatName=[], N_LEVELS=[], N_PROF=[], CYCLES=np.array([], dtype=int), floatType=[]
    )

    for v in ['PRES', 'TEMP', 'PSAL', 'SDN']:
        floatData[v] = np.array([])
        floatData[v + '_QC'] = np.array([])
    
    for v in ['WMO', 'LATITUDE', 'LONGITUDE', 'POSITION_QC', 'SDN_GRID', 'LATITUDE_GRID', 'LONGITUDE_GRID', 'CYCLE_GRID']:
        floatData[v] = np.array([])

    for v in ['DOXY', 'CHLA', 'BBP700', 'CDOM', 'NITRATE', 'DOWNWELLING_IRRADIANCE']:
        if v in common_variables:
            floatData[v] = np.array([])
            floatData[v + '_QC'] = np.array([])
            if v + '_ADJUSTED' in common_variables:
                floatData[v + '_ADJUSTED'] = np.array([])
                floatData[v + '_ADJUSTED' + '_QC'] = np.array([])

    for fn,cn in zip(files,core_files):
        print(fn, cn)
        # try to load the profile as absolute path or relative path
        try:
            nc = Dataset(fn, 'r')
        except:
            try:
                nc = Dataset(Path(ARGO_PATH) / fn, 'r')
            except:
                raise FileNotFoundError('No such file {} or {}'.format(fn, str(Path(ARGO_PATH) / fn)))

        try:
            cc = Dataset(cn, 'r')
        except:
            try:
                cc = Dataset(Path(ARGO_PATH) / cn, 'r')
            except:
                warnings.warn('No such file {} or {}'.format(fn, str(Path(ARGO_PATH) / fn)))

        # number of profile cycles
        M = cc.dimensions['N_LEVELS'].size
        N = cc.dimensions['N_PROF'].size

        wmo = ''
        if N > 1:
            for let in nc.variables['PLATFORM_NUMBER'][:][0,:].compressed():
                wmo = wmo + let.decode('UTF-8')
        else:
            for let in nc.variables['PLATFORM_NUMBER'][:].compressed():
                wmo = wmo + let.decode('UTF-8')

        cycle = nc.variables['CYCLE_NUMBER'][:].data.flatten()

        ftype = ''
        if 'PLATFORM_TYPE' in nc.variables.keys():
            for let in nc.variables['PLATFORM_TYPE'][:].compressed():
                ftype = ftype + let.decode('UTF-8')

        floatData['floatName']  = floatData['floatName'] + [int(wmo)]
        floatData['N_LEVELS']   = floatData['N_LEVELS']  + [M]
        floatData['N_PROF']     = floatData['N_PROF']    + [N]
        floatData['CYCLES']     = np.append(floatData['CYCLES'], cycle)
        floatData['CYCLE_GRID'] = np.append(floatData['CYCLE_GRID'], np.array(N*M*[cycle[0]]))
        floatData['floatType']  = floatData['floatType'] + [ftype]
        floatData['WMO']        = np.append(floatData['WMO'], np.array(M*N*[wmo]))

        # load in variables that will be in every file
        floatData['PRES'] = np.append(floatData['PRES'], cc.variables['PRES'][:].data.flatten())
        floatData['PRES_QC'] = np.append(floatData['PRES_QC'], read_qc(cc.variables['PRES_QC'][:].data.flatten()))
        floatData['TEMP'] = np.append(floatData['TEMP'], cc.variables['TEMP'][:].data.flatten())
        floatData['TEMP_QC'] = np.append(floatData['TEMP_QC'], read_qc(cc.variables['TEMP_QC'][:].data.flatten()))
        floatData['PSAL'] = np.append(floatData['PSAL'], cc.variables['PSAL'][:].data.flatten())
        floatData['PSAL_QC'] = np.append(floatData['PSAL_QC'], read_qc(cc.variables['PSAL_QC'][:].data.flatten()))
        floatData['SDN'] = np.append(floatData['SDN'], cc.variables['JULD'][:].data.flatten() + mdates.datestr2num('1950-01-01'))
        floatData['SDN_QC'] = np.append(floatData['SDN_QC'], read_qc(cc.variables['JULD_QC'][:].data.flatten()))
        floatData['SDN_GRID'] = np.append(floatData['SDN_GRID'], np.array(N*M*[np.nanmean(cc.variables['JULD'][:].data.flatten() + mdates.datestr2num('1950-01-01'))]))
        floatData['LATITUDE'] = np.append(floatData['LATITUDE'], cc.variables['LATITUDE'][:].data.flatten())
        floatData['LATITUDE_GRID'] = np.append(floatData['LATITUDE_GRID'], np.array(N*M*[np.nanmean(cc.variables['LATITUDE'][:].data.flatten())]))
        floatData['LONGITUDE'] = np.append(floatData['LONGITUDE'], cc.variables['LONGITUDE'][:].data.flatten())
        floatData['LONGITUDE_GRID'] = np.append(floatData['LONGITUDE_GRID'], np.array(N*M*[np.nanmean(cc.variables['LONGITUDE'][:].data.flatten())]))
        floatData['POSITION_QC'] = np.append(floatData['POSITION_QC'], read_qc(cc.variables['POSITION_QC'][:].data.flatten()))

        # loop through other possible BGC variables
        bgc_vars = ['DOXY', 'CHLA', 'BBP700', 'CDOM', 'NITRATE', 'DOWNWELLING_IRRADIANCE']
        for v in bgc_vars:
            if v in common_variables:
                floatData[v] = np.append(floatData[v], vertically_align(cc.variables['PRES'][:].data.flatten(), nc.variables['PRES'][:].data.flatten(), nc.variables[v][:].data.flatten()))
            v_adj = v + '_ADJUSTED'
            if v_adj in common_variables:
                floatData[v_adj] = np.append(floatData[v_adj], vertically_align(cc.variables['PRES'][:].data.flatten(), nc.variables['PRES'][:].data.flatten(), nc.variables[v_adj][:].data.flatten()))

        floatData['dPRES'] = delta_pres(cc.variables['PRES'][:].data.flatten(), nc.variables['PRES'][:].data.flatten())

        for v in floatData.keys():
            v_qc = v + '_QC'
            if v_qc in common_variables:
                floatData[v_qc] = np.append(floatData[v_qc], read_qc(nc.variables[v_qc][:].data.flatten()))

        if 'DOXY' in floatData.keys():
            floatData['O2Sat'] = 100*floatData['DOXY']/unit.oxy_sol(floatData['PSAL'], floatData['TEMP'], unit='micromole/kg')
            floatData['O2Sat_QC'] = get_worst_flag(floatData['TEMP_QC'], floatData['PSAL_QC'], floatData['DOXY_QC'])

    return floatData

def read_history_qctest(nc):

    QC_ACTION = np.squeeze(nc.variables['HISTORY_ACTION'][:].data)
    actions = []
    for row in QC_ACTION:
        rval = ''
        for let in row:
            rval = rval + let.decode('UTF-8')
        actions.append(rval.strip())
    actions = np.array(actions)

    QC_TESTS  = np.squeeze(nc.variables['HISTORY_QCTEST'][:].data)
    tests = []
    for row in QC_TESTS:
        rval = ''
        for let in row:
            rval = rval + let.decode('UTF-8')
        tests.append(rval.strip())
    tests = np.array(tests)

    qcp_index = np.logical_or(actions == 'QCP', actions == 'QCP$')
    qcf_index = np.logical_or(actions == 'QCF', actions == 'QCF$')
    QCP, QCF = tests[qcp_index][0], tests[qcf_index][0]

    return QCP, QCF

def dict_clean(float_data, bad_flags=None):

    clean_float_data = copy.deepcopy(float_data)
    qc_flags = [k for k in clean_float_data.keys() if '_QC' in k]

    if bad_flags is None:
        for qc_key in qc_flags:
            data_key   = qc_key.replace('_QC','')
            good_index = np.logical_or(np.logical_or(clean_float_data[qc_key] < 4, clean_float_data[qc_key] == 5), clean_float_data[qc_key] == 8)
            bad_index  = np.invert(good_index)

            if data_key == 'POSITION':
                for dk in ['LATITUDE', 'LONGITUDE']:
                    clean_float_data[dk][bad_index] = np.nan
            else:
                clean_float_data[data_key][bad_index] = np.nan
    else:
        if type(bad_flags) is int:
            bad_flags = [bad_flags]
        
        for flag in bad_flags:
            for qc_key in qc_flags:
                data_key = qc_key.replace('_QC','')
                bad_index = clean_float_data[qc_key] == flag

                if data_key == 'POSITION':
                    for dk in ['LATITUDE', 'LONGITUDE']:
                        clean_float_data[dk][bad_index] = np.nan
                else:
                    clean_float_data[data_key][bad_index] = np.nan
        
    return clean_float_data

def dict_fillvalue_clean(float_data):

    clean_float_data = copy.deepcopy(float_data)
    qc_keys = [k for k in clean_float_data.keys() if '_QC' in k and 'SDN' not in k]

    for k in qc_keys:
        data_key   = k.replace('_QC','')
        if data_key == 'POSITION':
            for dk in ['LATITUDE', 'LONGITUDE']:
                fillvalue_index = clean_float_data[dk] >= 99999. # use greater than because date fillval is 999999
            clean_float_data[dk][fillvalue_index] = np.nan
        else:
            fillvalue_index = clean_float_data[data_key] >= 99999. # use greater than because date fillval is 999999
            clean_float_data[data_key][fillvalue_index] = np.nan

    # check if there is in-air data present
    if 'PPOX_DOXY' in float_data.keys():
        fillvalue_index = clean_float_data['PPOX_DOXY'] >= 99999. # use greater than because date fillval is 999999
        clean_float_data['PPOX_DOXY'][fillvalue_index] = np.nan

    fillvalue_index = clean_float_data['SDN'] >= 999999.
    clean_float_data['SDN'][fillvalue_index] = np.nan

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


def calc_gain(data, ref, inair=True, zlim=25., verbose=True):
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

        intersect_cycles = np.intersect1d(cycle, np.unique(inair_cycle), assume_unique=True)

        mean_float_data = np.nan*np.ones((ref.shape[0],4))
        for i,c in enumerate(intersect_cycles):
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
        grid_cycle = data['CYCLE_GRID'][surf_ix]
        grid_time  = data['SDN_GRID'][surf_ix]
        cycle = data['CYCLES']
        time  = data['SDN']

        z_woa = ref['z']
        woa_data = ref['WOA']

        woa_index = np.where(z_woa <= zlim)[0]
        woa_surf = np.nanmean(woa_data[woa_index,:],axis=0)
        woa_surf = woa_data[0,:]

        mean_float_data = np.nan*np.ones((woa_surf.shape[0],4))
        g = np.nan*np.ones((woa_surf.shape[0],))
        for i,t in enumerate(time): # uncomment when ready
        # for i,c in enumerate(cycle):
            ref_o2sat = woa_surf[i]
            # subset_o2sat = surf_o2sat[grid_cycle == c]
            subset_o2sat = surf_o2sat[grid_time == t] # uncomment when ready
            mean_float_data[i,0] = c
            mean_float_data[i,1] = np.sum(~np.isnan(subset_o2sat))
            mean_float_data[i,2] = np.nanmean(subset_o2sat)
            mean_float_data[i,3] = np.nanstd(subset_o2sat)

            g[i] = ref_o2sat/mean_float_data[i,2]
        
        g[g == 0] = np.nan

        return g, mean_float_data, woa_surf

def grid_var(gridded_cycle, Nprof, Nlevel, argo_var):

    return gV

def vertically_align(P1, P2, V2):

	out = np.nan*np.ones(P1.shape)

	for i, p in enumerate(P1):
		index  = np.abs(P2 - p) == np.min(np.abs(P2 - p))
		out[i] = np.nanmean(V2[index])

	return out

def delta_pres(P1, P2):

	dpres = np.nan*np.ones(P1.shape)

	for i, p in enumerate(P1):
		index    = np.abs(P2 - p) == np.min(np.abs(P2 - p))
		dpres[i] = np.nanmean(P2[index] - p)

	return dpres

def doxy_range_check(floatdict, verbose=True, adjusted=False):
    cleandict = copy.deepcopy(floatdict)
    if adjusted:
        key = 'DOXY_ADJUSTED'
    else:
        key = 'DOXY'
    doxy = floatdict[key]
    outside_range = np.logical_or(doxy < -5, doxy > 600)
    if verbose:
        sys.stdout.write('{} valyes found outside RTQC range check, replacing with NaN'.format(np.sum(outside_range)))

    doxy[outside_range] = np.nan
    cleandict[key] = doxy

    return cleandict

def aic(data, resid):
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


def bic(data, resid):
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
    variable (v2) in a float_data dictionary (output of load_argo), though
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
    lut_data = np.loadtxt(REF_PATH / 'T_lL_tau_3830_4330.dat')
    lut_lL = lut_data[0,1:]
    lut_T  = lut_data[1:,0]
    tau100 = lut_data[1:,1:]
    thickness = thickness*np.ones((N-1,))

    # translate boundary layer thickness to temperature dependent tau
    f_thickness = RectBivariateSpline(lut_T, lut_lL, tau100, kx=1, ky=1)
    tau_T = np.squeeze(f_thickness(thickness, mean_temp, grid=False))

    # loop through oxygen data
    for i in range(N-1):
        dt = t_sec[i+1] - t_sec[i]

        # do the correction using the mean filter, get the mean time
        mean_oxy[i]  = (1/(2*oxy_b(dt, tau_T[i])))*(DO[i+1] - oxy_a(dt, tau_T[i])*DO[i])
    
    # interpolate back to original times for output
    f = interp1d(mean_time, mean_oxy, kind='linear', bounds_error=False, fill_value='extrapolate')
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
    f = interp1d(mean_time, mean_oxy, kind='linear', bounds_error=False, fill_value='extrapolate')
    DO_out = f(t_sec)

    return DO_out