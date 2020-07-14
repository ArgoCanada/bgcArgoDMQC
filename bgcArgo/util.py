#!/usr/bin/python

import numpy as np

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

        if lat_ix[-1] != lat_ix.shape[0] - 1:
            lat_ix = np.append(lat_ix, np.array([lat_ix[-1]+1]))

        return lat_ix

def get_lon_index(lon, lon_bounds, cross180):
    '''function to pull appropriate WOA longitude values, handles crossing 180'''
    
    if cross180:
        lon_ix = np.logical_or(lon <= lon_bounds[0], lon >= lon_bounds[1])
        lon_ix = np.where(lon_ix)[0]
        diff_index = np.where(np.diff(lon_ix) != 1)[0][0]
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
