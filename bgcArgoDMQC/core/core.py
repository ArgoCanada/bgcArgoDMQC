import sys
import copy
from pathlib import Path

import numpy as np
import matplotlib.dates as mdates

import gsw
from netCDF4 import Dataset

from .. import io
from .. import interp
from .. import unit

# ----------------------------------------------------------------------------
# LOCAL MACHINE SETUP
# ----------------------------------------------------------------------------

def set_dirs(**kwargs):
    '''
    Convenience wrapper for bgcArgoDMQC.io.PathHandler().set_dirs()

    Set local directories to look for Argo, WOA, and NCEP data.

    Args:
        argo_path (str or path-like): location of local Argo data
        ncep_data (str or path-like): location of local NCEP data
        woa_path (str or path-like): location of local World Ocean Atlas data
    '''

    io.Path.set_dirs(**kwargs)

def get_index(index='bgc', **kwargs):
    '''
    Get the global, biogeochemical, synthetic, or metadata Argo index. 

    Args:
        index (str): *bgc* for the biogeochemical Argo index, *global* for the core index, *synthetic* for the synthetic index, or *meta* for the metadata index
    '''
    if index == 'bgc':
        if '__bgcindex__' not in globals():
            global __bgcindex__
            __bgcindex__ = io.read_index()
        return_index = __bgcindex__
    elif index == 'global':
        if '__globalindex__' not in globals():
            global __globalindex__
            __globalindex__ = io.read_index(mission='C')
        return_index = __globalindex__
    elif index == 'synthetic':
        if '__synthindex__' not in globals():
            global __synthindex__
            __synthindex__ = io.read_index(mission='S')
        return_index = __synthindex__
    elif index == 'meta':
        if '__metaindex__' not in globals():
            global __metaindex__
            __metaindex__ = io.read_index(mission='M')
        return_index = __metaindex__
    elif index == 'traj':
        if '__trajindex__' not in globals():
            global __trajindex__
            __trajindex__ = io.read_index(mission='T')
        return_index = __trajindex__
    else:
        raise ValueError(f'Input "{index}" is unrecognized')

    for arg, val in kwargs.items():
        return_index = return_index[return_index[arg] == val]
    
    return return_index.reset_index()

# ----------------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------------

def load_argo(local_path, wmo, grid=False, verbose=True):
    '''
    Function to load in all data from a single float, using BRtraj, meta,
    and Sprof files.
    
    Args:
        local_path: local path of float data
        wmo: float ID number
    
    Returns:
        floatData: python dict() object with Argo variables
    
        CYCLES, LATITUDE, LONGITUDE, and SDN all also have
        analogous <VAR>_GRID fields that match the    
        dimension of PRES, TEMP, PSAL, DOXY, and O2SAT  
    
    Author:   
        Christopher Gordon
        Fisheries and Oceans Canada
        chris.gordon@dfo-mpo.gc.ca
    '''

    # make local_path a Path() object from a string, account for windows path
    local_path = Path(local_path)
    dac = io.get_dac(wmo)

    wmo = str(wmo) if type(wmo) is not str else wmo

    # check that necessary files exist - can continue without BRtraj file but
    # need Sprof and meta files
    BRtraj = local_path / dac / wmo / f'{wmo}_BRtraj.nc'
    Sprof  = local_path / dac / wmo / f'{wmo}_Sprof.nc'
    meta   = local_path / dac / wmo / f'{wmo}_meta.nc'

    # check if BRtraj is there, flag for moving forward if not
    BRtraj_flag = True
    if not BRtraj.exists():
        BRtraj_nc = None
        BRtraj_flag = False
        if verbose:
            sys.stdout.write('Continuing without BRtraj file\n')
    elif BRtraj.exists():
        BRtraj_nc = Dataset(BRtraj, 'r')
        if 'PPOX_DOXY' not in BRtraj_nc.variables.keys() and 'DOXY' not in BRtraj_nc.variables.keys():
            BRtraj_flag = False
            if verbose:
                sys.stdout.write('BRtraj file exists, but no in-air data exists, continuing without using BRtraj file\n')
    else:
        BRtraj_nc = None

    # Sprof and meta are required, so raise error if they are not there
    if not Sprof.exists():
        raise FileNotFoundError(f'No such Sprof file: {Sprof.absolute()}')
    if not meta.exists():
        raise FileNotFoundError(f'No such meta file: {meta.absolute()}')

    # load synthetic and meta profiles
    Sprof_nc = Dataset(Sprof, 'r')
    meta_nc  = Dataset(meta, 'r')

    # number of profile cycles
    M = Sprof_nc.dimensions['N_LEVELS'].size
    N = Sprof_nc.dimensions['N_PROF'].size

    # fillvalue dict
    fillvalue = {k:Sprof_nc[k]._FillValue for k in Sprof_nc.variables.keys()}
    
    floatData = read_flat_variables(Sprof_nc)
    floatData['SDN']  = floatData['JULD'] + mdates.datestr2num('1950-01-01')
    floatData['CYCLES'] = floatData['CYCLE_NUMBER']
    floatData['WMO'] = wmo

    qc_keys = [s for s in floatData.keys() if '_QC' in s and 'PROFILE' not in s]
    for qc in qc_keys:
        floatData[qc] = io.read_qc(floatData[qc])

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
        floatData['PDEN'] = gsw.pot_rho_t_exact(gsw.SA_from_SP(floatData['PSAL'], floatData['PRES'], floatData['LONGITUDE_GRID'], floatData['LATITUDE_GRID']), floatData['TEMP'], floatData['PRES'], 0)

    if 'DOXY' in floatData.keys():
        optode_flag = get_optode_type(int(wmo)) == 'AANDERAA_OPTODE_4330'
        floatData['O2Sat'] = 100*floatData['DOXY']/unit.oxy_sol(floatData['PSAL'], floatData['TEMP'], floatData['PDEN'], a4330=optode_flag)
        # match the fill values
        ix = np.logical_or(np.logical_or(floatData['PSAL'] == fillvalue['PSAL'], floatData['TEMP'] == fillvalue['TEMP']), floatData['DOXY'] == fillvalue['DOXY'])
        floatData['O2Sat'][ix] = fillvalue['DOXY']
        # get the worst QC flag from each quantity that goes into the calculation
        floatData['O2Sat_QC'] = copy.deepcopy(floatData['DOXY_QC'])

    if BRtraj_flag:
        if 'PPOX_DOXY' in BRtraj_nc.variables.keys() and 'TEMP_DOXY' in BRtraj_nc.variables.keys():
            floatData['PPOX_DOXY']  = BRtraj_nc.variables['PPOX_DOXY'][:].data.flatten()
            floatData['TEMP_DOXY']  = BRtraj_nc.variables['TEMP_DOXY'][:].data.flatten()
            floatData['TRAJ_CYCLE'] = BRtraj_nc.variables['CYCLE_NUMBER'][:].data.flatten()
            floatData['inair']      = True
        elif 'DOXY' in BRtraj_nc.variables.keys() and 'TEMP_DOXY' in BRtraj_nc.variables.keys():
            # unit conversion from umol kg-1 to pO2, some shaky S and P assumptions?
            floatData['PPOX_DOXY'] = unit.doxy_to_pO2(unit.umol_per_sw_to_mmol_per_L(
                BRtraj_nc.variables['DOXY'][:].data.flatten(),
                0, # salinity is 0 in air???
                BRtraj_nc.variables['TEMP_DOXY'][:].data.flatten(),
                0 # pressure is 0 in air???
            ), 0, BRtraj_nc.variables['TEMP_DOXY'][:].data.flatten())
            floatData['TEMP_DOXY']  = BRtraj_nc.variables['TEMP_DOXY'][:].data.flatten()
            floatData['TRAJ_CYCLE'] = BRtraj_nc.variables['CYCLE_NUMBER'][:].data.flatten()
            floatData['inair']      = True
        else:
            floatData['inair']      = False
    else:
        floatData['inair']          = False


    return floatData, Sprof, BRtraj, meta, fillvalue

def load_profile(local_path, wmo, cyc, kind='C', direction='A'):
    '''
    Function to load in all data from a single profile file,
    core or BGC.
    
    Args:
        local_path: local path of float data
        wmo: float ID number
        cyc: cycle number
        kind: core ("C") or B ("B") file
        direction: ascending ("A") or descending ("D")
    
    Returns:
        floatData: python dict() object with Argo variables
    
        CYCLES, LATITUDE, LONGITUDE, and SDN all also have
        analogous <VAR>_GRID fields that match the    
        dimension of PRES, TEMP, PSAL, DOXY, and O2SAT  
    
    Author:   
        Christopher Gordon
        Fisheries and Oceans Canada
        chris.gordon@dfo-mpo.gc.ca
    '''

    # make local_path a Path() object from a string, account for windows path
    local_path = Path(local_path)
    dac = io.get_dac(wmo)

    wmo = str(wmo) if type(wmo) is not str else wmo
    cyc = str(cyc) if type(wmo) is not str else cyc

    kind = '' if kind == 'C' else kind
    direction = '' if direction == 'A' else direction

    # check that the file exists - check for D-mode file first
    profFile = local_path / dac / wmo / 'profiles' / f'{kind}D{wmo}_{cyc:03d}{direction}.nc'
    profFile = profFile.parent / f'{kind}R{wmo}_{cyc:03d}{direction}.nc' if not profFile.exists() else profFile

    if not profFile.exists():
        raise FileNotFoundError(f'No R- or D-mode file: {profFile.absolute()}')
    
    nc = Dataset(profFile, 'r')

    # fillvalue dict
    fillvalue = {k:nc[k]._FillValue for k in nc.variables.keys()}
    
    floatData = read_flat_variables(nc)
    floatData['SDN']  = floatData['JULD'] + mdates.datestr2num('1950-01-01')
    floatData['CYCLES'] = floatData['CYCLE_NUMBER']
    floatData['WMO'] = wmo

    qc_keys = [s for s in floatData.keys() if '_QC' in s and ('PROFILE' not in s and 'HISTORY' not in s)]
    for qc in qc_keys:
        floatData[qc] = io.read_qc(floatData[qc])

    return floatData, profFile, fillvalue

def read_flat_variables(nc):
    '''
    Read all variables and dimensions from an Argo netCDF file.

    Args:
        nc: a netCDF file object
    
    Returns:
        floatData: python dict with all variable and dimension names
    '''

    floatData = dict()
    for name, dim in nc.dimensions.items():
        floatData[name] = dim.size

    n_prof = nc.dimensions['N_PROF'].size
    n_levels = nc.dimensions['N_LEVELS'].size
    n_prof_array = []
    n_levels_array = []
    for i in range(n_prof):
        n_prof_array = n_prof_array + n_levels*[i]
        n_levels_array = n_levels_array + list(range(n_levels))

    floatData['N_PROF_DIM'] = floatData['N_PROF']
    floatData['N_LEVELS_DIM'] = floatData['N_LEVELS']
    floatData['N_PROF'] = np.array(n_prof_array)
    floatData['N_LEVELS'] = np.array(n_levels_array)

    for name, var in nc.variables.items():
        floatData[name] = var[:].data.flatten()

    return floatData

def read_gridded_variables(nc):
    '''
    Read all variables and dimensions from an Argo Sprof file, do not flatten
    arrays, keep as 2D arrays.

    Args:
        nc: a netCDF file object
    
    Returns:
        floatData: python dict with all variable and dimension names
    '''

    floatData = dict()
    for name, dim in nc.dimensions.items():
        floatData[name] = dim.size
    for name, var in nc.variables.items():
        floatData[name] = var[:].data

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
    qc_flags = [s for s in clean_float_data.keys() if '_QC' in s and ('PROFILE' not in s and 'HISTORY' not in s)]

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
    qc_keys = [s for s in clean_float_data.keys() if '_QC' in s and ('PROFILE' not in s and 'HISTORY' not in s)]

    for k in qc_keys:
        data_key   = k.replace('_QC','')
        if data_key == 'POSITION':
            for dk in ['LATITUDE', 'LONGITUDE', 'LATITUDE_GRID', 'LONGITUDE_GRID']:
                if dk in clean_float_data.keys():
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

    if 'SDN_GRID' in float_data.keys():
        fillvalue_index = clean_float_data['SDN_GRID'] >= 999999.
        clean_float_data['SDN_GRID'][fillvalue_index] = np.nan

    return clean_float_data

def track(float_data):
    # make 'track' array with columns (time, lat, lon) to be used in interpolation
    track = np.array([float_data['SDN'], float_data['LATITUDE'], float_data['LONGITUDE']]).T

    return track

def woa_to_float_track(track, param, zlim=(0,1000), local_path='./', verbose=True):
    '''
    Function to load WOA18 climatological data for comparison with autonomous
    floats. Data to be interpolated along the provided track (t, lat, lon).
    Combines function load_woa_data() and interp_woa_data() for convenience,
    see documentation for those funcions for more detail.
    
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
        z: WOA depth array
        woa_interp: 2D array of requested WOA parameter (depth x time)
    
    Author:   
        Christopher Gordon
        Fisheries and Oceans Canada
        chris.gordon@dfo-mpo.gc.ca
    
    Last update: 2020-04-23
    
    Change log:
    '''

    xtrack, woa_track, woa_data = io.load_woa_data(track, param, zlim=zlim, local_path=local_path, verbose=verbose)
    woa_interp, wt, yrday = interp.interp_woa_data(xtrack, woa_track, woa_data, verbose=verbose)
    z = woa_track[0]

    return z, woa_interp, wt

def ncep_to_float_track(varname, track, local_path='./'):
    '''
    Function to load NCEP reanalysis data for comparison with autonomous
    floats. Data to be interpolated along the provided track (t, lat, lon).
    Combines function load_ncep_data() and interp_ncep_data() for convenience,
    see documentation for those funcions for more detail.
    
    Args:
        varname: either 'pres' (pressure) or 'rhum' (relative humidity)
        track: array with the columns (SDN, lat, lon)
    
    Returns:
        z: WOA depth array
        woa_interp: 2D array of requested WOA parameter (depth x time)
    
    Author:   
        Christopher Gordon
        Fisheries and Oceans Canada
        chris.gordon@dfo-mpo.gc.ca
    
    Last update: 2020-04-29
    
    Change log:
    '''

    xtrack, ncep_track, data = io.load_ncep_data(track, varname, local_path=local_path)
    if track[0,0] > ncep_track[0][-1] and mdates.num2date(track[0,0]).year == mdates.datetime.date.today().year:
        raise ValueError('First float date occurs after last NCEP date, NCEP data not available yet, recommend using WOA data to calcualte gain')
    ncep_interp, wt = interp.interp_ncep_data(xtrack, ncep_track, data)

    return ncep_interp, wt


def calc_gain(data, ref, inair=True, zlim=25., verbose=True):
    '''
    Calculate the gain for each profile by comparing float oxygen data to a
    reference data set, either NCEP for in-air or WOA surface data if in-air
    comparison is not available.
    
    Args:
        data: float data dict object, output from load_argo()
        ref: reference data set, either NCEP pO2 or WOA O2sat
        inair: boolean flag to indicate if comparison to NCEP in-air
            data or WOA surface data should be done, default to
            in-air, but function also performs check
        zlim: lower limit to define as 'surface' and take mean within,
            default value 25 dbar, for use only when inair is False
    
    Returns:
        g: vector of gains
        surf_data: array of float surface stats (cycle, N, mean, std)
    
    Author:   
        Christopher Gordon
        Fisheries and Oceans Canada
        chris.gordon@dfo-mpo.gc.ca
    
    Last update: 2020-04-23
    
    Change log:
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
        for i,t in enumerate(time):
            ref_o2sat = woa_surf[i]
            subset_o2sat = surf_o2sat[grid_time == t] # uncomment when ready
            mean_float_data[i,0] = cycle[i]
            mean_float_data[i,1] = np.sum(~np.isnan(subset_o2sat))
            mean_float_data[i,2] = np.nanmean(subset_o2sat)
            mean_float_data[i,3] = np.nanstd(subset_o2sat)

            g[i] = ref_o2sat/mean_float_data[i,2]
        
        g[g == 0] = np.nan

        return g, mean_float_data, woa_surf

def calc_gain_with_carryover(pO2_opt_air, pO2_ref_air, pO2_opt_water):
    '''
    Calculate gain with carryover parameter, following Bittig et al. (2018).

    Args:
        pO2_opt_air (array-like): partial pressure measured by the oxygen optode in-air
        pO2_ref_air (array-like): partial pressure in-air from a reference dataset such as NCEP
        pO2_opt_water (array-like): partial pressure of oxygen measured by the optode just below the surface

    Returns:
        *need to run this by Henry and see if I'm doing it right*

    Derive the O2 slope including a correction for 'carry-over' effect, to
    account for the observation that optode in-air data do not represent pure
    air but show a bias by in-water O2 saturation excess/deficiency (Bittig 
    and Kortzinger 2015). Johnson et al. (2015) confirm the 'carry-over' effect
    for optodes close to the surface (~20cm). 

    Carry-over effect is recommended to be account for Argo floats using in-air
    measurements, if enough surfacings are available (N > 20). It both removes
    an identified bias (which is most relevant for cases with strong 
    super-/undersaturation and/or carry-overs) and reduces uncertainty on the
    O2 slope factor. The equation for linear regression is as follows (see,
    e.g., Bittig et al., 2018):

    m*pO2^{optode}_{surf in-air} - pO2^{reference}_{in-air} 
        = c*(m*pO2^{optode}_{surf in-water} - pO2^{reference}_{in-air})

    where: 
        - m is the O2 slope factor: m = pO2_adjusted / pO2
        - pO2^{optode}_{surf in-air} is the oxygen partial pressure observed by
        the optode in-air (i.e., close to the water surface), e.g., MC = X+11
        - pO2^{reference}_{in-air} is the reference oxygen partial pressure in-air,
        e.g., from re-analysis data
        - pO2^{optode}_{surf in-water} is the oxygen partial pressure observed by 
        the optode at the water surface (in-water), e.g., MC = X+10 or profile 
        MC = X-10
        - c is the slope of the 'carry-over' effect, i.e., the water-fraction of 
        the observed optode in-air data.

    Above equation can be used for linear regression to obtain m and c from
    data of the partial pressures (from several cycles together). See 
    Thierry Virginie, Bittig Henry, The Argo-Bgc Team (2018). Argo quality 
    control manual for dissolved oxygen concentration. 
    https://doi.org/10.13155/46542
    '''

    x1 = pO2_opt_air - pO2_ref_air
    y1 = pO2_opt_water - pO2_ref_air

    x1 = x1[:,np.newaxis]

    carry_over, resid, _, _ = np.linalg.lstsq(x1, y1, rcond=None)
    c = carry_over

    gains = ((1-c)*pO2_ref_air)/(pO2_opt_air - c*pO2_opt_water)

    return gains, carry_over


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

def range_check(key, floatdict, verbose=True):
    if 'range_dict' not in globals():
        global range_dict
        range_dict = dict(
            PRES=(-5, np.inf),
            TEMP=(-2.5, 40),
            PSAL=(2, 41),
            DOXY=(-5, 600),
        )

    cleandict = copy.deepcopy(floatdict)

    argo_var = floatdict[key]
    r = range_dict[key.replace('_ADJUSTED','')]
    outside_range = np.logical_or(argo_var < r[0], argo_var > r[1])
    if verbose: # pragma: no cover
        sys.stdout.write('{} values found outside RTQC range check, replacing with NaN\n'.format(np.sum(outside_range)))

    argo_var[outside_range] = np.nan
    cleandict[key] = argo_var

    return cleandict

def calc_fixed_doxy_adjusted_error(S, T, P, fix_err=10):
    '''
    Calculate DOXY_ADJUSTED_ERROR for fixed partial pressure of 10 mbar 
    PPOX_DOXY.
    '''

    error = unit.pO2_to_doxy(np.array(S.shape[0]*[fix_err]), S, T, P=P)

    return error

def get_optode_type(wmo):
    if '__metaindex__' not in globals():
        global __metaindex__
        __metaindex__ = get_index(index='meta')
    
    ix = __metaindex__[__metaindex__.wmo == wmo]

    local_file = Path(io.Path.ARGO_PATH) / ix.dac.iloc[0] / str(wmo) / ix.file.iloc[0].split('/')[-1]
    nc = Dataset(local_file)

    doxy_index = io.get_parameter_index(nc['SENSOR'][:].data, 'OPTODE_DOXY')
    if doxy_index.shape[0] == 0:
        return 'NO_OPTODE_FOUND'
    else:
        optode_type = io.read_ncstr(nc['SENSOR_MODEL'][:].data[doxy_index[0], :])
        return optode_type
