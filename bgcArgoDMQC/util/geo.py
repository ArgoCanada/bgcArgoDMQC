
import numpy as np

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

def haversine(c1, c2):
    '''
    Calculate the distance between two coordinates using the Haversine formula.

    Args:
        c1 (tuple): first set of coordinates - (lat, long)
        c2 (tuple): second set of coordinates - (lat, long)

    Returns: 
        Distnace between c1 and c2 in kilometers
    '''

    # approx radius of earth in km
    R = 6373.0

    # convert coordinates to radians
    lat1 = np.radians(c1[0])
    lon1 = np.radians(c1[1])
    lat2 = np.radians(c2[0])
    lon2 = np.radians(c2[1])

    # difference in coords
    dlat = lat2 - lat1
    dlon = lon2 - lon1

    # Haversine formula
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    distance = R * c

    return distance
