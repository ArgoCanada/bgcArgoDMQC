
import numpy as np
import pandas as pd

from .core import *
from .. import unit
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

    THe main function serves to minimize the onus on the user to organize
    variables for quality control. Calculating an oxygen gain becomes simple::

        gains = syn.calc_gains(ref='NCEP')
    '''
    
    set_dirs = set_dirs

    def __init__(self, wmo, keep_fillvalue=False, rcheck=True, verbose=False):

        self.__floatdict__, self.__Sprof__, self.__BRtraj__, self.__meta__ = load_argo(ARGO_PATH, wmo, grid=True, verbose=verbose)
        self.__rawfloatdict__ = self.__floatdict__

        # local path info
        self.argo_path = ARGO_PATH
        self.woa_path  = WOA_PATH
        self.ncep_path = NCEP_PATH

        self.assign(self.__floatdict__)
        if not keep_fillvalue:
            self.rm_fillvalue()

        if rcheck:
            self.check_range('DOXY')

    def assign(self, floatdict):
        '''
        Assign variables from float dictionary (output of load_argo(...))
        to synthetic profile sprof() object.
        '''

        # metadata and dimension variables
        self.floatType  = floatdict['floatType']
        self.N_PROF     = floatdict['N_PROF']
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
        self.PDEN = gsw.pot_rho_t_exact(gsw.SA_from_SP(self.PSAL, self.PRES, self.LONGITUDE_GRID, self.LATITUDE_GRID), self.TEMP, self.PRES, 0) - 1000

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

    def get_gridded_var(self, *args):

        if not hasattr(self, '__griddict__'):
            self.__griddict__ = read_sprof_gridded_variables(Dataset())

        varlist = list()
        for v in args:
            varlist.append(self.__griddict__[v])

        return varlist

    def rm_fillvalue(self):
        '''
        Remove FillValue from all variables.
        '''
        self.__nofillvaluefloatdict__ = dict_fillvalue_clean(self.__rawfloatdict__)
        self.__floatdict__ = copy.deepcopy(self.__nofillvaluefloatdict__)
        self.assign(self.__nofillvaluefloatdict__)
        self.to_dataframe()

    def clean(self, bad_flags=None):
        '''
        Remove bad data from all variables, using <PARAM>_QC to determine bad data. 
        Optional input `bad_flags` can be used to specify which flag values are bad,
        with a default bad flags set to be 4, 6, 7.
        '''
        self.__cleanfloatdict__ = dict_clean(self.__floatdict__, bad_flags=bad_flags)
        self.__floatdict__ = copy.deepcopy(self.__cleanfloatdict__)
        self.assign(self.__cleanfloatdict__)
        self.to_dataframe()

    def reset(self):
        '''
        Reset all variables back to original loaded variables. Undoes the effect of
        clean(), rm_fillvalue(), check_range().
        '''
        self.__floatdict__ = copy.deepcopy(self.__rawfloatdict__)
        self.assign(self.__rawfloatdict__)
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

        self.assign(self.__rangecheckdict__)
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

        exvals = [k for k,v in self.__floatdict__.items() if type(v) is not np.ndarray]
        for key in list(set(self.__floatdict__.keys()) - set(df.columns) - set(exvals)):
            if self.__floatdict__[key].shape[0] == df.shape[0]:
                df[key] = self.__floatdict__[key]

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

        if not hasattr(self, 'track'):
            self.get_track()

        self.NCEP, self.__NCEPweights__ = ncep_to_float_track('pres', self.track, local_path=self.ncep_path)
        
        return copy.deepcopy(self.NCEP)

    def get_woa(self, verbose=True):
        '''
        Loads WOA data along the float track
        '''

        if not hasattr(self, 'track'):
            self.get_track()
        
        self.z_WOA, self.WOA, self.__WOAweights__ = woa_to_float_track(self.track, 'O2sat', local_path=self.woa_path, verbose=verbose)

        return copy.deepcopy(self.WOA)

    def calc_gains(self, ref='NCEP', zlim=25., verbose=True):
        '''
        Calculate gain values using NCEP or WOA reference data. Uses function
        calc_gain(). 
        '''

        if not hasattr(self, 'track'):
            self.get_track()

        if ref == 'NCEP':
            # check if reference data is already calculated
            if not hasattr(self, 'NCEP'):
                self.get_ncep(verbose=verbose)

            pH2O = unit.pH2O(util.get_var_by('TEMP_DOXY', 'TRAJ_CYCLE', self.__floatdict__))

            common_cycles, c1, c2 = np.intersect1d(self.CYCLE, np.unique(self.__floatdict__['TRAJ_CYCLE']), assume_unique=True, return_indices=True)

            self.NCEP_PPOX = unit.atmos_pO2(self.NCEP[c1], pH2O[c2])/100
            self.__NCEPgains__, self.__NCEPfloatref__ = calc_gain(self.__floatdict__, self.NCEP_PPOX, verbose=verbose)
            self.gains = self.__NCEPgains__

        if ref == 'WOA':
            # check if reference data is already calculated
            if not hasattr(self, 'WOA'):
                self.get_woa(verbose=verbose)

            self.__WOAgains__, self.__WOAfloatref__, self.__WOAref__ = calc_gain(self.__floatdict__, dict(z=self.z_WOA, WOA=self.WOA), inair=False, zlim=zlim, verbose=verbose)
            self.gains = self.__WOAgains__
        
        self.gain = np.nanmean(self.gains)
        
        return copy.deepcopy(self.gains)

    def calc_fixed_error(self, fix_err=10):

        self.DOXY_ADJUSTED_ERROR = calc_fixed_doxy_adjusted_error(self.__floatdict__, fix_err=fix_err)
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

            if not hasattr(self, 'df'):
                self.to_dataframe()

            g = plot.var_cscatter(self.df, varname=var, **kwargs)

        elif kind == 'profiles':
            varlist = kwargs.pop('varlist')

            if not hasattr(self, 'df'):
                self.to_dataframe()

            g = plot.profiles(self.df, varlist=varlist, **kwargs)

        elif kind == 'qcprofiles':
            varlist = kwargs.pop('varlist')

            if not hasattr(self, 'df'):
                self.to_dataframe()

            g = plot.qc_profiles(self.df, varlist=varlist, **kwargs)

        else:
            raise ValueError('Invalid input for keyword argument "kind"')

        return g


    def describe(self):
        '''
        Describe the dataframe of data stored in the sprof object.
        '''

        if not hasattr(self, 'df'):
            self.to_dataframe()
        
        print('Data for synthetic profile file for float {}'.format(self.WMO))

        sys.stdout.write('Variables:\n')
        for k in self.__floatdict__.keys():
            sys.stdout.write('{}\n'.format(k))
        sys.stdout.write('\n')
    
    def update_field(self, field, value, where=None):

        where = slice(None) if where is None else where
        self.__floatdict__[field][where] = value

        if field == 'DOXY' or field == 'DOXY_QC':
            optode_flag = get_optode_type(int(self.__floatdict__['WMO'])) == 'AANDERAA_OPTODE_4330'
            self.__floatdict__['O2Sat'] = unit.oxy_saturation(self.__floatdict__['DOXY'], self.__floatdict__['PSAL'], self.__floatdict__['TEMP'], self.__floatdict__['PDEN'], a4330=optode_flag)
            self.__floatdict__['O2Sat_QC'] = util.get_worst_flag(self.__floatdict__['TEMP_QC'], self.__floatdict__['PSAL_QC'], self.__floatdict__['DOXY_QC'])

        self.assign(self.__floatdict__)
        self.to_dataframe()
    
    def export_files(self, data_mode='D', glob=None):

        glob = 'BR*.nc' if glob is None else glob
        r_files = (self.__Sprof__.parent / 'profiles').glob(glob)

        io.export_files(self.__floatdict__, r_files, self.gain, data_mode=data_mode)

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

        if type(date) is str:
            date = mdates.datestr2num(date)
        if type(date) is mdates.datetime.datetime:
            date = mdates.date2num(date)

        meta_dict = dict(date=date)
        if lat is None:
            lat = np.nan
        if lon is None:
            lon = np.nan
        meta_dict['lat'] = lat
        meta_dict['lon'] = lon

        if data_dict is not None and len(kwargs) > 0:
            # apppend kwargs to dict
            for k in kwargs.keys():
                data_dict[k] = kwargs.pop(k)
        data_dict = dict(**kwargs) if data_dict is None else data_dict
        

        # default label value        
        if label is None:
            label = ' '
        # if there isn't already independent data, make a dict for it
        if not hasattr(self, '__indepdict__'):
            self.__indepdict__ = {label:data_dict}
            self.__indepmeta__ = {label:meta_dict}
        # if there is one already, append to it
        else:
            self.__indepdict__[label] = data_dict
            self.__indepmeta__[label] = meta_dict

    def compare_independent_data(self, fmt='*', **kwargs):
        '''
        Plot the independent data overtop of the nearest profile in time
        '''

        if not hasattr(self, 'df'):
            self.df = self.to_dataframe()

        plot_dict = copy.deepcopy(self.__indepdict__)
        meta_dict = copy.deepcopy(self.__indepmeta__)

        fig, ax_list = plot.compare_independent_data(self.df, plot_dict, meta_dict, fmt, **kwargs)
        return fig, ax_list