
import numpy as np
import pandas as pd

from .core import *
from .. import unit
from .. import util
from .. import plot
from .. import io

class sprof:
    '''
    Class that loads Argo synthetic profile data for a given float ID number
    (wmo). 

    Then, load the individual variables into fields in the class, for
    example::

        syn = sprof(wmo)
        print(syn.DOXY)

    Or load it into a pandas dataframe::

        df = syn.to_dataframe()

    The main function serves to minimize the onus on the user to organize
    variables for quality control. Calculating an oxygen gain becomes simple::

        gains = syn.calc_gains(ref='NCEP')
    '''
    
    def __init__(self, wmo, keep_fillvalue=False, rcheck=True, verbose=False):

        self.__floatdict__, self.__Sprof__, self.__BRtraj__, self.__meta__, self.__fillvalue__ = load_argo(io.Path.ARGO_PATH, wmo, grid=True, verbose=verbose)
        self.__rawfloatdict__ = self.__floatdict__

        # local path info
        self.argo_path = io.Path.ARGO_PATH
        self.woa_path  = io.Path.WOA_PATH
        self.ncep_path = io.Path.NCEP_PATH

        self.WMO = wmo

        self.get_track()
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
        priority_vars = ['PRES', 'PRES_QC', 'TEMP', 'TEMP_QC', 'PSAL', 'PSAL_QC']
        bgc_vars = list(set(self.__floatdict__.keys()) & set(['DOXY', 'DOXY_QC', 'DOXY_ADJUSTED', 'DOXY_ADJUSTED_QC', 'CHLA', 'CHLA_QC', 'CHLA_ADJUSTED', 'CHLA_ADJUSTED_QC', 'BBP700', 'BBP700_QC', 'BBP700_ADJUSTED', 'BBP_ADJUSTED_QC']))
        priority_vars = priority_vars + bgc_vars

        for v in priority_vars:
            df[v] = self.__floatdict__[v]

        n_level = df.shape[0]

        for k in set(self.__floatdict__.keys()) - set(df.columns):
            dim = self.__floatdict__[k].shape[0] if type(self.__floatdict__[k]) is np.ndarray else np.inf
            if dim == n_level:
                df[k.replace('_GRID', '')] = self.__floatdict__[k]
        
        self.df = df

        return copy.deepcopy(self.df)

    def get_track(self):
        '''
        Creates a track array with columns::

            [serial datenum, latitude, longitude]

        the track array is used for the interpolation of reference data along
        the float track.
        '''
        self.track = track(self.__floatdict__)

        return copy.deepcopy(self.track)

    def get_ncep(self, verbose=True):
        '''
        Loads NCEP data along the float track
        '''

        self.NCEP, self.__NCEPweights__ = ncep_to_float_track('pres', self.track, local_path=self.ncep_path)
        
        return copy.deepcopy(self.NCEP)

    def get_woa(self, verbose=True):
        '''
        Loads WOA data along the float track
        '''

        self.z_WOA, self.WOA, self.__WOAweights__ = woa_to_float_track(self.track, 'O2sat', local_path=self.woa_path, verbose=verbose)

        return copy.deepcopy(self.WOA)

    def calc_gains(self, ref='NCEP', zlim=25., verbose=True):
        '''
        Calculate gain values using NCEP or WOA reference data. Uses function
        calc_gain(). 
        '''

        if ref == 'NCEP':
            pH2O = unit.pH2O(util.get_var_by('TEMP_DOXY', 'TRAJ_CYCLE', self.__floatdict__))
            _, c1, c2 = np.intersect1d(self.CYCLE_NUMBER, np.unique(self.__floatdict__['TRAJ_CYCLE']), assume_unique=True, return_indices=True)

            try:
                self.NCEP_PPOX = unit.atmos_pO2(self.NCEP[c1], pH2O[c2])/100
            except KeyError:
                self.get_ncep(verbose=verbose)
                self.NCEP_PPOX = unit.atmos_pO2(self.NCEP[c1], pH2O[c2])/100

            self.__NCEPgains__, self.__NCEPfloatref__ = calc_gain(self.__floatdict__, self.NCEP_PPOX, verbose=verbose)
            self.gains = self.__NCEPgains__

        if ref == 'WOA':
            try:
                self.__WOAgains__, self.__WOAfloatref__, self.__WOAref__ = calc_gain(self.__floatdict__, dict(z=self.z_WOA, WOA=self.WOA), inair=False, zlim=zlim, verbose=verbose)
            except KeyError:
                self.get_woa()
                self.__WOAgains__, self.__WOAfloatref__, self.__WOAref__ = calc_gain(self.__floatdict__, dict(z=self.z_WOA, WOA=self.WOA), inair=False, zlim=zlim, verbose=verbose)
            self.gains = self.__WOAgains__
        
        self.gain = np.nanmean(self.gains)
        
        return copy.deepcopy(self.gains)

    def calc_fixed_error(self, fix_err=6):

        self.DOXY_ADJUSTED_ERROR = calc_fixed_doxy_adjusted_error(self.PSAL, self.TEMP, self.PRES, fix_err=fix_err)
        self.__floatdict__['DOXY_ADJUSTED_ERROR'] = self.DOXY_ADJUSTED_ERROR

        return copy.deepcopy(self.DOXY_ADJUSTED_ERROR)

    def plot(self, kind, **kwargs):

        if kind == 'gain':
            ref = kwargs['ref']
            if ref == 'NCEP':
                g = plot.gainplot(self.SDN, self.__NCEPfloatref__[:,2], self.NCEP_PPOX, self.__NCEPgains__, ref)
            elif ref == 'WOA':
                g = plot.gainplot(self.SDN, self.__WOAfloatref__[:,2], self.__WOAref__, self.__WOAgains__, ref)
            else:
                raise ValueError('Invalid input for keyword argument "ref"')

        elif kind == 'cscatter':
            var = kwargs.pop('varname')

            g = plot.variable_color_scatter(self.df, varname=var, **kwargs)

        elif kind == 'profiles':
            varlist = kwargs.pop('varlist')

            g = plot.profiles(self.df, varlist=varlist, **kwargs)

        elif kind == 'qcprofiles':
            varlist = kwargs.pop('varlist')

            g = plot.qc_profiles(self.df, varlist=varlist, **kwargs)

        else:
            raise ValueError('Invalid input for keyword argument "kind"')

        return g


    def describe(self):
        '''
        Describe the dataframe of data stored in the sprof object.
        '''

        return self.df.describe()
    
    def update_field(self, field, value, where=None):

        current_float_dict = copy.deepcopy(self.__floatdict__)
        self.reset()

        where = slice(None) if where is None else where
        self.__floatdict__[field][where] = value

        if field in ['DOXY', 'TEMP', 'PSAL']:
            optode_flag = get_optode_type(int(self.__floatdict__['WMO'])) == 'AANDERAA_OPTODE_4330'
            self.__floatdict__['O2Sat'] = unit.oxy_saturation(self.__floatdict__['DOXY'], self.__floatdict__['PSAL'], self.__floatdict__['TEMP'], self.__floatdict__['PDEN'], a4330=optode_flag)
        elif field == 'DOXY_QC':
            self.__floatdict__['O2Sat_QC'] = copy.deepcopy(self.__floatdict__['DOXY_QC'])
        
        self.__floatdict__ = current_float_dict
        self.to_dataframe()

    def set_fillvalue(self, field, where=None):

        self.update_field(field, self.__fillvalue__[field], where)
    
    def export_files(self, data_mode='D', glob=None, **kwargs):

        current_float_dict = copy.deepcopy(self.__floatdict__)
        self.reset()

        glob = 'BR*.nc' if glob is None else glob
        files = (self.__Sprof__.parent / 'profiles').glob(glob)
        self.reset()

        io.export_files(self.__floatdict__, files, self.gain, data_mode=data_mode, **kwargs)
        self.__floatdict__ = current_float_dict

    def add_independent_data(self, date=None, lat=None, lon=None, data_dict=None, label=None, **kwargs):
        '''
        Add independent data in order to easily plot and compare.

        Args:
            date (optional, str, float, or date):
                Date as a string ('YYYY-MM-DD'), serial datenumber, or python
                date object. 
            data_dict(optional, dict): 
                Dictionary containing variables to be added, key names should
                match the naming convention to Argo variables
            label (optional, str):
                Label or name of the dataset being added, important when adding
                more than one source of independent data
            **kwargs:
                Argo variable names and their values, essentailly the same as
                inputting a dictionary without having to actually build one

        Returns:
            None
        
        Example::
            syn = sprof(wmo_number)
            df = pd.read_csv('my_winkler_data.csv')

            syn.add_independent_data(
                date='2020-10-04', # metadata arguments, optional, if no date matches to first profile
                label='Winkler' # label to classify the data - for if you have more than one source
                PRES=df.pressure, # data arguments, match naming to Argo variables
                DOXY=df.dissolved_oxygen,
                LATITUDE=df.lat, LONGITUDE=df.lon, # again, optional
            )

            data = dict(PRES=df.PRES, DOXY=df.DOXY)
            syn.add_independent_data(data_dict=data, date='2020-10-04')
        '''

        
        date = mdates.datestr2num(date) if type(date) is str else date
        date = mdates.date2num(date) if type(date) is mdates.datetime.datetime else date
        meta_dict = dict(date=date)

        lat = np.nan if lat is None else lat
        lon = np.nan if lon is None else lon
        meta_dict['lat'] = lat
        meta_dict['lon'] = lon

        if data_dict is not None and len(kwargs) > 0:
            raise ValueError('Cannot input data as both dict and kwargs')
        data_dict = dict(**kwargs) if data_dict is None else data_dict
        

        # default label value        
        if label is None:
            label = ' '
            
        try:
            self.__indepdict__[label] = data_dict
            self.__indepmeta__[label] = meta_dict
        except KeyError:
            self.__indepdict__ = {label:data_dict}
            self.__indepmeta__ = {label:meta_dict}

    def compare_independent_data(self, fmt='*'):
        '''
        Plot the independent data overtop of the nearest profile in time
        '''

        plot_dict = copy.deepcopy(self.__indepdict__)
        meta_dict = copy.deepcopy(self.__indepmeta__)

        print(self.df, self.WMO, plot_dict, meta_dict, fmt)

        g = plot.compare_independent_data(self.df, self.WMO, plot_dict, meta_dict, fmt=fmt)
        return g