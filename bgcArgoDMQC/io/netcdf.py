#!/usr/bin/python
import sys

from pathlib import Path
import numpy as np
import pandas as pd
from netCDF4 import Dataset

from ..util import refill_array
from ..configure import read_config

def read_ncstr(arr):
    decode_str = np.array([f.decode('utf-8') for f in arr])
    out = ''.join(decode_str)

    return out.strip()

def read_qc(flags):

    decode_flags = np.array([f.decode('utf-8') for f in flags])
    decode_flags[decode_flags == ' '] = '4'

    out_flags = np.array([int(f) for f in decode_flags])

    return out_flags

def get_parameter_index(parameter_array, parameter):
    str_arr = np.array([read_ncstr(s) for s in parameter_array])
    index = np.where(str_arr == parameter)[0]

    return index

def create_fillvalue_array(nc_var):

    new_var = np.full(
        nc_var.shape,
        nc_var._FillValue,    
        dtype=nc_var.datatype
    )

    return new_var

def string_to_array(s, dim, encode='utf-8'):
    '''
    Args:
        s: input string or comment
        dim: netCDF dimension which the string will be made to fill using trailing whitespace
        encode (optional): encoding each character, default utf-8 (def. python netCDF format)

    Returns:
        numpy.array of single letters

    Author:
        Christopher Gordon
        Fisheries and Oceans Canada
        chris.gordon@dfo-mpo.gc.ca

    Change log:
        - 2021-07-05: initial commit
    '''

    M = dim.size - len(s)
    full_string = s + M*' '
    str_array = np.array([f'{c}'.encode(encode) for c in full_string])

    return str_array

def copy_netcdf(infile, outfile, exclude_vars=[], exclude_dims=[]):
    '''
    Copy data from a netCDF file with the exception of dimension and variable
    names listed in exclude_vars and exclude_dims.
    '''
    with Dataset(infile) as src, Dataset(outfile, 'w') as dst:
        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)
        # copy dimensions except for the excluded
        for name, dimension in src.dimensions.items():
            if name not in exclude_dims:
                dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
        # copy file data except for the excluded
        for name, variable in src.variables.items():
            if name not in exclude_vars:
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                # copy variable attributes all at once via dictionary
                dst[name].setncatts(src[name].__dict__)
                dst[name][:] = src[name][:]
    
    return Dataset(outfile, 'r+')

def iterate_dimension(infile, outfile, iterated_dimension, n=1):
    '''
    Add 1 or optional argument `n` to the size limited dimension `dimension`
    and fill the new dimension with the proper FillValue for each variable.
    '''

    with Dataset(infile) as src, Dataset(outfile, 'w') as dst:
        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)
        # copy dimensions except for the excluded
        for name, dimension in src.dimensions.items():
            if iterated_dimension == name:
                dst.createDimension(name, (len(dimension)+n if not dimension.isunlimited() else None))
            else:
                dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
        for name, variable in src.variables.items():
            new_variable = dst.createVariable(name, variable.datatype, variable.dimensions)
            # copy variable attributes all at once via dictionary
            dst[name].setncatts(src[name].__dict__)
            if iterated_dimension in variable.dimensions:
                # add FillValues along the added dimesion size
                arr = create_fillvalue_array(new_variable)
                # get location of iterated dimension
                axis = variable.dimensions.index(iterated_dimension)
                n_index = len(variable.dimensions)
                dst[name][:] = refill_array(axis, n_index, arr, src[name][:])
            else:
                # fill the variable with the same data as before
                dst[name][:] = src[name][:]

    return Dataset(outfile, 'r+')

def append_variable(fn, *args):
    '''
    Add an arbitrary number of variables (*args) to the existing netcdf file
    input fn. The input structure for each variable should be a dictionary with
    fields that can be passed directly to the netCDF file. If the variable
    name already exists, it will overwrite it with the new information.

    Args: 
        fn: string pointing to netcdf (.nc) file to be appended
        *args: arbitrary number of python dicts with all required fields to create
        or overwrite a new netcdf variable. Example: 

        new_var = dict(
            name='MY_NEW_VARIABLE',     # variable name, can be new or existing
            datatype=np.float64,        # variable datatype, can be np datatype or string
            dimensions=('N', 'M'),      # will be created with warning if they do not exist
            data=data_arr,              # the data
            long_name='The new variable',
            standard_name='my_new_var',
            units='degree_celsius',
            valid_min=0,
            valid_max=1e9,
            resolution=0.001,
            comment='Added by John Doe on Sept 20, 2019'
        )
    '''

    nc = Dataset(fn, 'r+')
    
    for new_var in args:
        name = new_var.pop('name')
        data = new_var.pop('data')
        nc.createVariable(name, new_var.pop('datatype'), new_var.pop('dimensions'))
        nc[name].setncatts(new_var)
        nc[name][:] = data

    return nc

def update_history(nc, dct):
    '''
    Update HISTORY_<PARAM> values in an Argo netCDF file
    '''

    hix = nc.dimensions['N_HISTORY'].size
    for name, value in dct.items():
        for i in range(nc.dimensions['N_PROF'].size):
            nc[name][hix,i,:] = string_to_array(value, nc.dimensions[nc[name].dimensions[-1]])

def profile_qc(flags):
    '''
    Return overall profile quality flag via the following from the Argo User
    Manual (v 3.41):

    3.2.2 Reference table 2a: overall profile quality flag
    https://vocab.nerc.ac.uk/collection/RP2/current
    N is defined as the percentage of levels with good data where:
    - QC flag values of 1, 2, 5, or 8 are considered GOOD data
    - QC flag values of 9 (missing) or “ “ are NOT USED in the computation
    All other QC flag values are BAD data
    The computation should be taken from <PARAM_ADJUSTED>_QC if available and from 
    <PARAM>_QC otherwise.
    n Meaning
    "" No QC performed
    A N = 100%; All profile levels contain good data.
    B 75% <= N < 100%
    C 50% <= N < 75%
    D 25% <= N < 50%
    E 0% < N < 25%
    F N = 0%; No profile levels have good data.

    Args:
        - flags (pandas.Series): quality flags for a given profile
    Returns:
        - grade (str): profile grade based on description above
    '''
    
    n_good = flags.isin([1, 2, 5, 8]).sum()
    n_exclude = flags.isin([9]).sum()

    pct = 100*n_good/(flags.size - n_exclude)

    grade = np.nan

    if flags.isin([0]).sum() >= flags.size - n_exclude:
        grade = ''

    if pct == 100:
        grade = 'A'
    elif pct >= 75:
        grade = 'B'
    elif pct >= 50:
        grade = 'C'
    elif pct >= 25:
        grade = 'D'
    elif pct > 0:
        grade = 'E'
    elif pct == 0:
        grade = 'F'

    if not type(grade) == str and np.isnan(grade):
        raise ValueError('No grade assigned, check input value of `flags`')

    return grade

def export_files(fdict, r_files, gain, data_mode='D', comment=None, equation=None, coeff=None):

    config = read_config()
    dmqc_date = pd.Timestamp.now(tz='utc').strftime('%Y%m%d%H%M%S')

    for fn in r_files:
        # define path to file, make directory if it does not exist
        dac = fn.as_posix().split('/')[-4]
        D_file = Path(fn.as_posix().replace('BR', f'B{data_mode}').\
            replace(f'dac/{dac}/', f'dac/{dac}/D/'))
        if not D_file.parent.exists():
            D_file.parent.mkdir(parents=True)
        sys.stdout.write(f'Working on D-mode file {D_file.as_posix()}...')

        D_nc = copy_netcdf(fn, D_file)
        last_calib = D_nc.dimensions['N_CALIB'].size-1

        # index for this cycle
        cycle = int(fn.as_posix().split('_')[-1].split('.')[0].replace('D', ''))
        ix = fdict['CYCLE_GRID'] == cycle
        N = D_nc.dimensions['N_LEVELS'].size

        # find index for DOXY along PARAMETER
        doxy_index = [read_ncstr(a) for a in D_nc['PARAMETER'][:].data[0,0,:,:]].index('DOXY')

        # fill in string info
        temp_comment  = f'Oxygen gain calculated following Johnson et al. 2015, doi:10.1175/JTECH-D-15-0101.1, using comparison between float and WOA data. Adjustment applied by {config["operator"]} ({config["affiliation"]}, orcid: {config["orcid"]})'
        comment  = temp_comment if comment is None else comment
        equation = 'DOXY_ADJUSTED = G*DOXY' if equation is None else equation
        coeff    = f'G = {gain:f}' if coeff is None else coeff
 
        # apply info to all profiles in file (not sure if this would ever not apply 
        # take caution when N_PROF > 1)
        for i in range(D_nc.dimensions['N_PROF'].size):
            D_nc['SCIENTIFIC_CALIB_COMMENT'][i,last_calib,doxy_index,:] = string_to_array(comment, D_nc.dimensions['STRING256'])
            D_nc['SCIENTIFIC_CALIB_EQUATION'][i,last_calib,doxy_index,:] = string_to_array(equation, D_nc.dimensions['STRING256'])
            D_nc['SCIENTIFIC_CALIB_COEFFICIENT'][i,last_calib,doxy_index,:] = string_to_array(coeff, D_nc.dimensions['STRING256'])

        D_nc['DOXY_QC'][:] = fdict['DOXY_QC'][ix][:N]
        D_nc['DOXY_ADJUSTED'][:] = fdict['DOXY_ADJUSTED'][ix][:N]
        D_nc['DOXY_ADJUSTED_QC'][:] = fdict['DOXY_ADJUSTED_QC'][ix][:N]
        D_nc['DOXY_ADJUSTED_ERROR'][:] = fdict['DOXY_ADJUSTED_ERROR'][ix][:N]

        for i in range(D_nc.dimensions['N_PROF'].size):
            flags = read_qc(D_nc['DOXY_ADJUSTED_QC'][:].data[i,:])
            grade = profile_qc(pd.Series(flags)).encode('utf-8')
            D_nc['PROFILE_DOXY_QC'][i] = grade
        
        data_state_indicator = create_fillvalue_array(D_nc['DATA_STATE_INDICATOR'])
        for i in range(D_nc.dimensions['N_PROF'].size):
            data_state_indicator[i,:] = string_to_array('2C+', D_nc.dimensions['STRING4'])
        D_nc['DATA_STATE_INDICATOR'][:] = data_state_indicator

        nc_data_mode = create_fillvalue_array(D_nc['DATA_MODE'])
        for i in range(D_nc.dimensions['N_PROF'].size):
            nc_data_mode[i] = data_mode
        D_nc['DATA_MODE'][:] = nc_data_mode

        parameter_data_mode = create_fillvalue_array(D_nc['PARAMETER_DATA_MODE'])
        for i in range(D_nc.dimensions['N_PROF'].size):
            tmp_pdm = D_nc['PARAMETER_DATA_MODE'][:].data[i,:]
            tmp_pdm[get_parameter_index(D_nc['PARAMETER'][:][i,0,:,:].data, 'DOXY')] = data_mode
            parameter_data_mode[i,:] = tmp_pdm
        D_nc['PARAMETER_DATA_MODE'][:] = parameter_data_mode


        history_dict = dict(
            HISTORY_INSTITUTION='BI',
            HISTORY_STEP='ARSQ',
            HISTORY_SOFTWARE='BGQC',
            HISTORY_SOFTWARE_RELEASE='v0.2',
            HISTORY_DATE=dmqc_date,
            HISTORY_ACTION='O2QC'
        )
        D_nc['DATE_UPDATE'][:] = string_to_array(dmqc_date, D_nc.dimensions['DATE_TIME'])

        update_history(D_nc, history_dict)
        sys.stdout.write('done\n')
        D_nc.close()