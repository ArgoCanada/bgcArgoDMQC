
import numpy as np
import pandas as pd

from .core import *
from .. import unit
from .. import util
from .. import plot
from .. import io

class bio_prof:
    '''
    Class that loads Argo synthetic profile data for a given float ID number
    (wmo). 

    Then, load the individual variables into fields in the class, for
    example::

        profiles = bio_prof(wmo)
        # dataframe-like, can export to actual df using profiles.to_dataframe()
        print(profiles.DOXY)
    '''
    
    set_dirs = set_dirs

    def __init__(self, wmo, cycles=None, keep_fillvalue=False, rcheck=True, verbose=False):

        self.__floatdict__ = {}

        # self.__floatdict__, self.__Sprof__, self.__BRtraj__, self.__meta__, self.__fillvalue__ = load_argo(ARGO_PATH, wmo, grid=True, verbose=verbose)
        self.__rawfloatdict__ = self.__floatdict__

        # local path info
        self.argo_path = io.Path.ARGO_PATH
        self.woa_path  = io.Path.WOA_PATH
        self.ncep_path = io.Path.NCEP_PATH

        # self.to_dataframe()

        if not keep_fillvalue:
            self.rm_fillvalue()

        if rcheck:
            self.check_range('DOXY')

    def __getitem__(self, index):
        return pd.Series(self.__floatdict__[index])

    def __setitem__(self, index, value):
        self.df[index] = value

    def __getattr__(self, index):
        return pd.Series(self.__floatdict__[index])

    def rm_fillvalue(self):
        '''
        Remove FillValue from all variables.
        '''
        self.__nofillvaluefloatdict__ = dict_fillvalue_clean(self.__rawfloatdict__)
        self.__floatdict__ = copy.deepcopy(self.__nofillvaluefloatdict__)
        self.to_dataframe()
    
    def clean(self, bad_flags=None):
        '''
        Remove bad data from all variables, using <PARAM>_QC to determine bad data. 
        Optional input `bad_flags` can be used to specify which flag values are bad,
        with a default bad flags set to be 4, 6, 7.
        '''
        self.__cleanfloatdict__ = dict_clean(self.__floatdict__, bad_flags=bad_flags)
        self.__floatdict__ = copy.deepcopy(self.__cleanfloatdict__)
        self.to_dataframe()

    def reset(self):
        '''
        Reset all variables back to original loaded variables. Undoes the effect of
        clean(), rm_fillvalue(), check_range().
        '''
        self.__floatdict__ = copy.deepcopy(self.__rawfloatdict__)
        self.to_dataframe()

    def check_range(self, key, verbose=False):
        '''
        Performs a range check for variables that have a RTQC range available.
        Replaces values outside the range with NaN values. Takes string input
        to do the range check on that variable. Available variables are
        PRES, TEMP, PSAL, and DOXY. Can also take input 'all' to do the range
        check on all four variables, or a list of strings to do each of those
        variables.
        '''
        if key == 'all':
            key = ['PRES', 'TEMP', 'PSAL', 'DOXY']
        elif type(key) is not list:
            key = [key]
        
        for k in key:
            if k in self.__floatdict__.keys():
                self.__rangecheckdict__ = range_check(k, self.__floatdict__, verbose=verbose)
                self.__floatdict__ = self.__rangecheckdict__

                # recalculate O2sat if its DOXY
                if k == 'DOXY':
                    optode_flag = get_optode_type(int(self.__rangecheckdict__['WMO'])) == 'AANDERAA_OPTODE_4330'
                    self.__rangecheckdict__['O2Sat'] = 100*self.__rangecheckdict__['DOXY']/unit.oxy_sol(self.__rangecheckdict__['PSAL'], self.__rangecheckdict__['TEMP'], self.__rangecheckdict__['PDEN'], a4330=optode_flag)

        self.to_dataframe()
    
    def to_dict(self):
        '''
        Returns a deepcopy of __floatdict__, which is the currect active
        dictionary (i.e. subject to the effects of clean(), reset(), etc.)
        '''
        return copy.deepcopy(self.__floatdict__)
    
    def to_dataframe(self):
        '''
        Returns a pandas dataframe containing data from the synthetic oxygen
        profile file.
        '''

        df = pd.DataFrame()
        n_level = self.__floatdict__['N_LEVELS']        
        priority_vars = ['PRES', 'PRES_QC', 'TEMP', 'TEMP_QC', 'PSAL', 'PSAL_QC']
        bgc_vars = list(set(self.__floatdict__.keys()) & set(['DOXY', 'DOXY_QC', 'DOXY_ADJUSTED', 'DOXY_ADJUSTED_QC', 'CHLA', 'CHLA_QC', 'CHLA_ADJUSTED', 'CHLA_ADJUSTED_QC', 'BBP700', 'BBP700_QC', 'BBP700_ADJUSTED', 'BBP_ADJUSTED_QC']))
        priority_vars = priority_vars + bgc_vars

        for v in priority_vars:
            df[v] = self.__floatdict__[v]

        for k in set(self.__floatdict__.keys()) - set(df.columns):
            dim = self.__floatdict__[k].shape[0] if type(self.__floatdict__[k]) is np.ndarray else np.inf
            if dim == n_level:
                df[k] = self.__floatdict__[k]

        return copy.deepcopy(self.df)
    
    def update_field(self, field, value, where=None):

        where = slice(None) if where is None else where
        self.__floatdict__[field][where] = value

        if field in ['DOXY', 'TEMP', 'PSAL']:
            optode_flag = get_optode_type(int(self.__floatdict__['WMO'])) == 'AANDERAA_OPTODE_4330'
            self.__floatdict__['O2Sat'] = unit.oxy_saturation(self.__floatdict__['DOXY'], self.__floatdict__['PSAL'], self.__floatdict__['TEMP'], self.__floatdict__['PDEN'], a4330=optode_flag)
        elif field == 'DOXY_QC':
            self.__floatdict__['O2Sat_QC'] = copy.deepcopy(self.__floatdict__['DOXY_QC'])

        self.assign(self.__floatdict__)
        self.to_dataframe()
    
    def set_fillvalue(self, field, where=None):

        self.update_field(field, self.__fillvalue__[field], where)
    
    def export_files(self, data_mode=None, glob=None):

        io.export_files(self.__floatdict__, self.__files__, data_mode=data_mode)