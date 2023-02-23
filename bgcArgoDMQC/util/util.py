import sys
from pathlib import Path
import numpy as np
from netCDF4 import Dataset

def decode_woa_var(varname):
    '''take WOA variable name input and output relevant info'''

    input_to_woa_param = dict(
        T='t',
        S='s',
        O2='o',
        O2sat='O',
        NO3='n',
        Si='i',
        PO4='p',
    )

    input_to_woa_ftype = dict(
        T='A5B7',
        S='A5B7',
        O2='all',
        O2sat='all',
        NO3='all',
        Si='all',
        PO4='all',
    )

    input_to_woa_folder = dict(
        T='temperature',
        S='salinity',
        O2='oxygen',
        O2sat='o2sat',
        NO3='nitrate',
        Si='silicate',
        PO4='phosphate',
    )

    param  = input_to_woa_param[varname]
    ftype  = input_to_woa_ftype[varname]
    ftpdir = input_to_woa_folder[varname]

    return param, ftype, ftpdir

def get_qctests(hex_code):

    # hex to numeric
    num = int(hex_code, 16)
    # list to save test number in
    tests = []
    for i in range(26,0,-1):
        qc_binary_id = 2**i
        if qc_binary_id <= num:
            num -= qc_binary_id
            tests.append(i)
    
    if num != 0:
        sys.stdout.write('NOTE: decoding QC tests left a non-zero remainder, suggest investigation\n')

    return tests[::-1]

def display_qctests(QCP, QCF):

    QCP_numbers = get_qctests(QCP)
    QCF_numbers = get_qctests(QCF)

    test_descriptions = [
        'Platform Identification test\t\t\t ',
        'Impossible Date test\t\t\t\t ',
        'Impossible Location test\t\t\t ',
        'Position on Land test\t\t\t\t ',
        'Impossible Speed test\t\t\t\t ',
        'Global Range test\t\t\t\t ',
        'Regional Global Parameter test\t\t ',
        'Pressure Increasing test\t\t\t ',
        'Spike test\t\t\t\t\t ',
        'Top and Bottom Spike test (obsolete)\t\t ',
        'Gradient test\t\t\t\t\t ',
        'Digit Rollover test\t\t\t\t ',
        'Stuck Value test\t\t\t\t ',
        'Density Inversion test\t\t\t ',
        'Grey List test\t\t\t\t ',
        'Gross Salinity or Temperature Sensor Drift test',
        'Visual QC test\t\t\t\t ',
        'Frozen profile test\t\t\t\t ',
        'Deepest pressure test\t\t\t\t ',
        'Questionable Argos position test\t\t ',
        'Near-surface unpumped CTD salinity test\t ',
        'Near-surface mixed air/water test\t\t '
    ]

    sys.stdout.write('\n---------------------------------------------------------------------------\n')
    sys.stdout.write('| Test\t| Pass/Fail\t| Test name\t\t\t\t\t  |\n')
    sys.stdout.write('---------------------------------------------------------------------------\n')

    for i,t in enumerate(test_descriptions):

        if i+1 in QCF_numbers:
            pfn = 'Failed'
            c = '\x1b[0;0;41m'
        elif i+1 in QCP_numbers:
            pfn = 'Passed'
            c = ''
        else:
            pfn = 'Not performed'
            c = '\x1b[0;30;43m'

        sys.stdout.write(f'{c}| {i+1:d}\t| {pfn}\t| {t} |\x1b[0m\n')

    sys.stdout.write('---------------------------------------------------------------------------\n')

def utf_decode(nc_arr, verbose=True):

    dlist = []
    for row in nc_arr:
        rval = ''
        for let in row:
            rval = rval + let.decode('UTF-8')
        
        if verbose:
            print(rval)

        dlist.append(rval.strip())

    return dlist

def read_gain_value(fn, verbose=True):

    nc = Dataset(fn, 'r')

    eq    = nc.variables['SCIENTIFIC_CALIB_EQUATION']
    coeff = nc.variables['SCIENTIFIC_CALIB_COEFFICIENT']
    comm  = nc.variables['SCIENTIFIC_CALIB_COMMENT']

    nprof =  nc.dimensions['N_PROF'].size
    ncalib =  nc.dimensions['N_CALIB'].size

    if nprof == 1 and ncalib == 1:
        eqs    = np.array(utf_decode(np.squeeze(eq[:].data), verbose=verbose)).flatten()
        coeffs = np.array(utf_decode(np.squeeze(coeff[:].data), verbose=verbose)).flatten()
        comms  = np.array(utf_decode(np.squeeze(comm[:].data), verbose=verbose)).flatten()

        ix = np.array(['DOXY_ADJUSTED' in s for s in eqs])
        if np.sum(ix) == 0:
            return np.array(['G=NaN']), np.array(['DOXY_ADJUSTED=DOXY*G']), np.array(['No gain value found'])
        else:
            G = coeffs[ix]
            equation = eqs[ix]
            comment = comms[ix]

    elif nprof > 1 and ncalib == 1:
        eqs    = np.array([utf_decode(np.squeeze(eqq), verbose=verbose) for eqq in eq[:].data]).flatten()
        coeffs = np.array([utf_decode(np.squeeze(cqq), verbose=verbose) for cqq in coeff[:].data]).flatten()
        comms  = np.array([utf_decode(np.squeeze(mqq), verbose=verbose) for mqq in comm[:].data]).flatten()

        ix = np.array(['DOXY_ADJUSTED' in s for s in eqs])
        if np.sum(ix) == 0:
            return np.array(['G=NaN']), np.array(['DOXY_ADJUSTED=DOXY*G']), np.array(['No gain value found'])
        else:
            G = coeffs[ix]
            equation = eqs[ix]
            comment = comms[ix]
        
    elif nprof == 1 and ncalib > 1:
        eqs    = np.array([utf_decode(np.squeeze(eqq), verbose=verbose) for eqq in np.squeeze(eq[:].data)]).flatten()
        coeffs = np.array([utf_decode(np.squeeze(cqq), verbose=verbose) for cqq in np.squeeze(coeff[:].data)]).flatten()
        comms  = np.array([utf_decode(np.squeeze(mqq), verbose=verbose) for mqq in np.squeeze(comm[:].data)]).flatten()

        ix = np.array(['DOXY_ADJUSTED' in s for s in eqs])
        if np.sum(ix) == 0:
            return np.array(['G=NaN']), np.array(['DOXY_ADJUSTED=DOXY*G']), np.array(['No gain value found'])
        else:
            G = coeffs[ix]
            equation = eqs[ix]
            comment = comms[ix]

    elif nprof > 1 and ncalib > 1:
        eqs    = np.array([utf_decode(np.squeeze(eqq), verbose=verbose) for p1 in eq[:].data for eqq in p1]).flatten()
        coeffs = np.array([utf_decode(np.squeeze(cqq), verbose=verbose) for p1 in coeff[:].data for cqq in p1]).flatten()
        comms  = np.array([utf_decode(np.squeeze(mqq), verbose=verbose) for p1 in comm[:].data for mqq in p1]).flatten()

        ix = np.array(['DOXY_ADJUSTED' in s for s in eqs])
        if np.sum(ix) == 0:
            return np.array(['G=NaN']), np.array(['DOXY_ADJUSTED=DOXY*G']), np.array(['No gain value found'])
        else:
            G = coeffs[ix]
            equation = eqs[ix]
            comment = comms[ix]
    
    else:
        raise ValueError('Cannot parse calibration info based on dimensions N_PROF={} and N_CALIB={}'.format(nprof, ncalib))
    

    return G, equation, comment

def get_var_by(v1, v2, float_data):
    '''
    Function to calculate the mean of one variable (v1) indexed by a second
    variable (v2) in a float_data dictionary (output of load_argo), though
    it would work with any python dict
    
    Args:
        v1 (str): string input of a key in float_data
        v2 (str): string input of a key in float_data
        float_data (dict): python dict() object
    
    Returns:
        1D numpy array with mean values
    '''

    
    index = np.unique(float_data[v2])
    out_array = np.nan*np.ones((len(index)))
    for i,v in enumerate(index):
        out_array[i] = np.nanmean(float_data[v1][float_data[v2] == v])
    
    return out_array

def dict_append(d1, d2, mismatch_ok=False):
    '''
    Append two dicts or the intersection of their keys. Values for each key
    must be compatible with `np.append`.
    '''

    if mismatch_ok:
        keys = list(set(d1.keys()).intersection(set(d2.keys())))
    else:
        keys = d1.keys()
    
    d_out = {}
    for k in keys:
        d_out[k] = np.append(d1[k], d2[k])

def cycle_from_time(time, t_arr, c_arr):
    ix = np.abs(t_arr - time) == np.nanmin(np.abs(t_arr - time))
    return c_arr.loc[ix].iloc[0]

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

def get_worst_flag(*args):
    out_flags = np.ones(args[0].shape)

    if len(args) == 1: # pragma: no cover
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
                elif f in [5,8] and out_flags[i] <= 2:
                    out_flags[i] = f
                elif f == 9:
                    out_flags[i] = 9
        
    return out_flags
