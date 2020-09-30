#!/usr/bin/python

import sys
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


def get_lat_index(lat, lat_bounds):
        '''function to pull appropriate WOA latitude values'''

        lat_ix = np.logical_and(lat >= lat_bounds[0], lat <= lat_bounds[1])

        if not lat_ix.any():
            lat_ix = np.where(np.abs(lat - np.mean(lat_bounds)) == np.min(np.abs(lat - np.mean(lat_bounds))))[0]
        else:
            lat_ix = np.where(lat_ix)[0]
            
        if lat_ix[0] != 0:
            lat_ix = np.append(np.array([lat_ix[0]-1]), lat_ix)

        if lat_ix[-1] != lat.shape[0] - 1:
            lat_ix = np.append(lat_ix, np.array([lat_ix[-1]+1]))

        return lat_ix

def get_lon_index(lon, lon_bounds, cross180):
    '''function to pull appropriate WOA longitude values, handles crossing 180'''
    
    if cross180:
        lon_ix = np.logical_or(lon <= lon_bounds[0], lon >= lon_bounds[1])
        lon_ix = np.where(lon_ix)[0]
        diff_index = np.where(np.diff(lon_ix) != 1)[0]
        if diff_index.shape[0] != 0:
            diff_index = diff_index[0]
            half1_lon_ix = np.append(lon_ix[:diff_index], np.array([lon_ix[diff_index]+1]))
            half2_lon_ix = np.append(np.array([lon_ix[diff_index+1] - 1]), lon_ix[diff_index+1:])
            lon_ix = np.append(half1_lon_ix, half2_lon_ix)

        if lon_ix[0] != 0:
            lon_ix = np.append(np.array([lon_ix[0]-1]), lon_ix)

        if lon_ix[-1] != lon.shape[0] - 1:
            lon_ix = np.append(lon_ix, np.array([lon_ix[-1]+1]))

    else:
        lon_ix = np.logical_and(lon >= lon_bounds[0], lon <= lon_bounds[1])
        if not lon_ix.any():
            lon_ix = np.where(np.abs(lon - np.mean(lon_bounds)) == np.min(np.abs(lon - np.mean(lon_bounds))))[0]
        else:
            lon_ix = np.where(lon_ix)[0]
            
        if lon_ix[0] != 0:
            lon_ix = np.append(np.array([lon_ix[0]-1]), lon_ix)

        if lon_ix[-1] != lon.shape[0] - 1:
            lon_ix = np.append(lon_ix, np.array([lon_ix[-1]+1]))

    return lon_ix

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

    sys.stdout.write('---------------------------------------------------------------------------\n')
    sys.stdout.write('| Test\t| Pass/Fail\t| Test name\t\t\t\t\t  |\n')
    sys.stdout.write('---------------------------------------------------------------------------\n')

    for i,t in enumerate(test_descriptions):

        if i+1 in QCF_numbers:
            pfn = 'Failed'
        elif i+1 in QCP_numbers:
            pfn = 'Passed'
        else:
            pfn = 'Not performed'

        sys.stdout.write('| {:d}\t| {}\t| {} |\n'.format(i+1, pfn, t))

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

def read_gain_value(fn):

    nc = Dataset(fn, 'r')

    eq    = nc.variables['SCIENTIFIC_CALIB_EQUATION']
    coeff = nc.variables['SCIENTIFIC_CALIB_COEFFICIENT']
    comm  = nc.variables['SCIENTIFIC_CALIB_COMMENT']

    nprof =  nc.dimensions['N_PROF'].size
    ncalib =  nc.dimensions['N_CALIB'].size

    if nprof == 1 and ncalib == 1:
        eqs    = np.array(utf_decode(np.squeeze(eq[:].data), verbose=True)).flatten()
        coeffs = np.array(utf_decode(np.squeeze(coeff[:].data))).flatten()
        comms  = np.array(utf_decode(np.squeeze(comm[:].data))).flatten()

        ix = np.array(['DOXY_ADJUSTED' in s for s in eqs])
        if np.sum(ix) == 0:
            return np.array([np.nan]), np.array(['No gain value found'])
        else:
            G = coeffs[ix]
            comment = comms[ix]

    elif nprof > 1 and ncalib == 1:
        eqs    = np.array([utf_decode(np.squeeze(eqq), verbose=True) for eqq in eq[:].data]).flatten()
        coeffs = np.array([utf_decode(np.squeeze(cqq)) for cqq in coeff[:].data]).flatten()
        comms  = np.array([utf_decode(np.squeeze(mqq)) for mqq in comm[:].data]).flatten()

        ix = np.array(['DOXY_ADJUSTED' in s for s in eqs])
        if np.sum(ix) == 0:
            return np.array([np.nan]), np.array(['No gain value found'])
        else:
            G = coeffs[ix]
            comment = comms[ix]
        
    elif nprof == 1 and ncalib > 1:
        eqs    = np.array([utf_decode(np.squeeze(eqq), verbose=True) for eqq in np.squeeze(eq[:].data)]).flatten()
        coeffs = np.array([utf_decode(np.squeeze(cqq)) for cqq in np.squeeze(coeff[:].data)]).flatten()
        comms  = np.array([utf_decode(np.squeeze(mqq)) for mqq in np.squeeze(comm[:].data)]).flatten()

        ix = np.array(['DOXY_ADJUSTED' in s for s in eqs])
        if np.sum(ix) == 0:
            return np.array([np.nan]), np.array(['No gain value found'])
        else:
            G = coeffs[ix]
            comment = comms[ix]

    elif nprof > 1 and ncalib > 1:
        eqs    = np.array([utf_decode(np.squeeze(eqq), verbose=True) for p1 in eq[:].data for eqq in p1]).flatten()
        coeffs = np.array([utf_decode(np.squeeze(cqq)) for p1 in coeff[:].data for cqq in p1]).flatten()
        comms  = np.array([utf_decode(np.squeeze(mqq)) for p1 in comm[:].data for mqq in p1]).flatten()

        ix = np.array(['DOXY_ADJUSTED' in s for s in eqs])
        if np.sum(ix) == 0:
            return np.array([np.nan]), np.array(['No gain value found'])
        else:
            G = coeffs[ix]
            comment = comms[ix]
    

    return G, comment