
import numpy as np
import pandas as pd

from .core import *
from .. import unit
from .. import util
from .. import plot
from .. import io

class prof:
    '''
    Class that loads Argo synthetic profile data for a given float ID number
    (wmo). 

    Then, load the individual variables into fields in the class, for
    example::

        p = prof(wmo, cyc)
        print(p.DOXY)

    Or load it into a pandas dataframe::

        df = p.to_dataframe()
    '''
    
    set_dirs = set_dirs

    def __init__(self, wmo, cycle, kind='C', keep_fillvalue=False, rcheck=True, verbose=False):

        self.__floatdict__, self.__prof__, self.__fillvalue__ = load_profile(io.Path.ARGO_PATH, wmo, cycle, kind=kind)
        self.__rawfloatdict__ = self.__floatdict__
        self._dict = 'raw'

        # local path info
        self.argo_path = io.Path.ARGO_PATH
        self.woa_path  = io.Path.WOA_PATH
        self.ncep_path = io.Path.NCEP_PATH

        self.WMO = wmo
        self.cycle = cycle
        self.kind = kind

        self.to_dataframe()

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
        self._dict = 'nofill'
        self.to_dataframe()
    
    def clean(self, bad_flags=None):
        '''
        Remove bad data from all variables, using <PARAM>_QC to determine bad data. 
        Optional input `bad_flags` can be used to specify which flag values are bad,
        with a default bad flags set to be 4, 6, 7.
        '''
        self.__cleanfloatdict__ = dict_clean(self.__floatdict__, bad_flags=bad_flags)
        self.__floatdict__ = copy.deepcopy(self.__cleanfloatdict__)
        self._dict = 'clean'
        self.to_dataframe()

    def reset(self):
        '''
        Reset all variables back to original loaded variables. Undoes the effect of
        clean(), rm_fillvalue(), check_range().
        '''
        self.__floatdict__ = copy.deepcopy(self.__rawfloatdict__)
        self._dict = 'raw'
        self.to_dataframe()
    
    def set_dict(self, name):
        '''
        Call on the proper dict function given a dictionary name
        '''

        if name == 'raw':
            self.reset()
        elif name == 'nofill':
            self.rm_fillvalue()
        elif name == 'clean':
            self.clean()
        else:
            raise ValueError('Unregognized dict type name, should be one of "clean", "nofill", "raw"')

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
                if k == 'DOXY' and self.kind == 'S':
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
        priority_vars = ['PRES', 'PRES_QC', 'TEMP', 'TEMP_QC', 'PSAL', 'PSAL_QC'] if self.kind in ['C', 'S'] else ['PRES']
        bgc_vars = list(set(self.__floatdict__.keys()) & set(['DOXY', 'DOXY_QC', 'DOXY_ADJUSTED', 'DOXY_ADJUSTED_QC', 'CHLA', 'CHLA_QC', 'CHLA_ADJUSTED', 'CHLA_ADJUSTED_QC', 'BBP700', 'BBP700_QC', 'BBP700_ADJUSTED', 'BBP_ADJUSTED_QC']))
        priority_vars = priority_vars + bgc_vars

        for v in priority_vars:
            df[v] = self.__floatdict__[v]

        for k in set(self.__floatdict__.keys()) - set(df.columns):
            dim = self.__floatdict__[k].shape[0] if type(self.__floatdict__[k]) is np.ndarray else np.inf
            if dim == n_level:
                df[k] = self.__floatdict__[k]

        self.df = df
        return copy.deepcopy(self.df)
    
    def update_field(self, field, value, where=None):

        current_float_dict = copy.deepcopy(self._dict)
        self.reset()

        where = slice(None) if where is None else where
        self.__floatdict__[field][where] = value

        self.set_dict(current_float_dict)
        self.to_dataframe()
    
    def set_fillvalue(self, field, where=None):

        self.update_field(field, self.__fillvalue__[field], where)
    
    def export_files(self):

        io.update_nc(self.__floatdict__, self.__prof__)