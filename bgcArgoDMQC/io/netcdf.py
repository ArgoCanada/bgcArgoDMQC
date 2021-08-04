#!/usr/bin/python

from pathlib import Path
from netCDF4 import Dataset

def populate_flags():
    '''
    Args:
        arg: description

    Returns:
        out: description

    Author:
        Christopher Gordon
        Fisheries and Oceans Canada
        chris.gordon@dfo-mpo.gc.ca

    Change log:
        - 2021-07-05: initial commit
    '''

def generate_comments_equations(variable, gain=None, operator='[operator name]', affiliation='[operator affiliation]', orcid=''):
    '''
    Args:
        arg: description

    Returns:
        out: description

    Author:
        Christopher Gordon
        Fisheries and Oceans Canada
        chris.gordon@dfo-mpo.gc.ca

    Change log:
        - 2021-07-05: initial commit
    '''

def adjust_flags(index, flags):
    '''
    Args:
        arg: description

    Returns:
        out: description

    Author:
        Christopher Gordon
        Fisheries and Oceans Canada
        chris.gordon@dfo-mpo.gc.ca

    Change log:
        - 2021-07-05: initial commit
    '''

def fill_string_array(s):
    '''
    Args:
        arg: description

    Returns:
        out: description

    Author:
        Christopher Gordon
        Fisheries and Oceans Canada
        chris.gordon@dfo-mpo.gc.ca

    Change log:
        - 2021-07-05: initial commit
    '''

### older function

def copy_netcdf_except(infile, outfile, exclude_vars=[], exclude_dims=[]):
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
    
    return Dataset(outfile, 'a')
        
def append_variable_to_file(fn, *args):
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

    nc = Dataset(fn, 'a')
    
    for new_var in args:
        name = new_var.pop('name')
        data = new_var.pop('data')
        nc.createVariable(name, new_var.pop('datatype'), new_var.pop('dimensions'))
        nc[name].setncatts(new_var)
        nc[name][:] = data

    return nc