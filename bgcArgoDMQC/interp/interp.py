import sys

import numpy as np
import matplotlib.dates as mdates
from scipy.interpolate import interp1d

def interp_ncep_data(track, ncep_track, data):
    '''
    Function to interpolate NCEP reanalysios data along the provided 
    track (t, lat, lon).
    
    INPUT:
              track: array with the columns (SDN, lat, lon)
              ncep_track: tuple with NCEP time, lat, lon arrays
              data: output array from load_ncep_data()
    
    OUTPUT:
              ncep_interp: 2D array of requested WOA parameter (depth x time)
        
    AUTHOR:   Christopher Gordon
              Fisheries and Oceans Canada
              chris.gordon@dfo-mpo.gc.ca
    
    ACKNOWLEDGEMENT: this code is adapted from the SOCCOM SAGE_O2Argo matlab
    code, available via https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING,
    written by Tanya Maurer & Josh Plant
    
    LAST UPDATE: 23-04-2020
    
    CHANGE LOG:
    '''
    # extract ncep variables
    ncep_time, lat_sub, lon_sub = ncep_track

    # get float lats and lons, round out times at endpoints
    xt   = np.ones((track.shape[0]+2))
    xlat = np.ones((track.shape[0]+2))
    xlon = np.ones((track.shape[0]+2))

    xt[1:-1]   = track[:,0]
    xlat[1:-1] = track[:,1]
    xlon[1:-1] = track[:,2]

    xt[0]    = np.floor(xt[1])
    xt[-1]   = np.ceil(xt[-2])
    xlat[0]  = xlat[1]
    xlat[-1] = xlat[-2]
    xlon[0]  = xlon[1]
    xlon[-1] = xlon[-2]

    # interpolate float to ncep time grid
    f = interp1d(xt, xlat, bounds_error=False)
    ilat = f(ncep_time)
    f = interp1d(xt, xlon, bounds_error=False)
    ilon = f(ncep_time)

    # remove NaN values
    ix = ~np.isnan(ilat)
    ilat = ilat[ix]
    ilon = ilon[ix]
    time_sub = ncep_time[ix]

    xwt = [[],[]]

    # loop through NCEP time and get data at each location
    ncep_data = np.nan*np.ones((time_sub.shape[0],))
    for i in range(time_sub.shape[0]):
        # temporal weighted average
        time_pt = np.where(ncep_time >= time_sub[i])[0][0]
        tmp = data[time_pt,:,:]
        
        # latitudinal weighted average
        lat_pt = np.where(lat_sub >= ilat[i])[0][-1]
        wt = (lat_sub[lat_pt] - ilat[i]) / (lat_sub[lat_pt] - lat_sub[lat_pt+1])
        tmp = wt*tmp[lat_pt+1,:] + (1-wt)*tmp[lat_pt,:]

        xwt[0].append(wt)

        # longitudinal weighted average - final pt. 
        lon_pt = np.where(lon_sub >= ilon[i])[0][0]
        wt = (lon_sub[lon_pt] - ilon[i]) / (lon_sub[lon_pt] - lon_sub[lon_pt-1])
        tmp = wt*tmp[lon_pt-1] + (1-wt)*tmp[lon_pt]

        xwt[1].append(wt)

        ncep_data[i] = tmp

    # interpolate back to float time
    f = interp1d(time_sub, ncep_data, bounds_error=False)
    ncep_interp = f(track[:,0])

    return ncep_interp, xwt

def interp_woa_data(track, woa_track, data, verbose=True):
    '''
    Function to interpolate WOA18 climatological data along the provided 
    track (t, lat, lon).
    
    INPUT:
              track: array with the columns (SDN, lat, lon)
              woa_track: tuple with WOA time, lat, lon arrays
              data: output array from load_woa_data()
    
    OUTPUT:
              woa_interp: 2D array of requested WOA parameter (depth x time)
        
    AUTHOR:   Christopher Gordon
              Fisheries and Oceans Canada
              chris.gordon@dfo-mpo.gc.ca
    
    ACKNOWLEDGEMENT: this code is adapted from the SOCCOM SAGE_O2Argo matlab
    code, available via https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING,
    written by Tanya Maurer & Josh Plant
    
    LAST UPDATE: 23-04-2020
    
    CHANGE LOG:
    '''

    # get relevant WOA information
    z, lat, lon = woa_track

    # put float and WOA on common time axis/year
    M = z.shape[0]
    N = track.shape[0]
    yrday = np.array([dn - mdates.datestr2num('{}-01-01'.format(mdates.num2date(dn).year)) for dn in track[:,0]])
    woa_yrday = np.array([mdates.datestr2num('2020-{:02d}-15'.format(i+1)) - mdates.datestr2num('2020-01-01') for i in range(12)])

    # array for output
    woa_interp = np.nan*np.ones((M,N))

    xwt = [[],[],[]]

    for i in range(N):
        # leave values as nan if lat/lon/time are nan:
        if any(np.isnan(track[i,:])):
            continue

        # ---------------------------------------------------------------------
        # Get bounding indices, interp weights, and subset before interpolation
        # ---------------------------------------------------------------------
        t1 = np.where(woa_yrday < yrday[i])[0]
        if t1.shape[0] == 0: # early Jan
            t1 = 0
            t2 = 11
            wt = (yrday[i] + 15) / 30
        elif t1[-1] == 11: # late Dec
            t1 = t1[-1]
            t2 = 0
            wt = (365 - yrday[i] + 15) / 30
        else:
            t1 = t1[-1]
            t2 = t1 + 1
            dx1 = woa_yrday[t2] - woa_yrday[t1]
            dx2 = yrday[i] - woa_yrday[t1]
            wt = (dx1 - dx2) / dx1

        xwt[0].append(wt)

        lat_ix1 = np.where(lat < track[i,1])[0][-1]
        lat_ix2 = lat_ix1 + 1
        if lat_ix2 == lat.shape[0]:
            if verbose:
                sys.stdout.write('NOTE: latitude is unbounded, giving all averaging weight to nearest observation\n')
            lat_ix1 -= 1
            lat_ix2 -= 1
            lat_wt = 0
        else:
            dx1 = lat[lat_ix2] - lat[lat_ix1]
            dx2 = track[i,1] - lat[lat_ix1]
            lat_wt = (dx1 - dx2) / dx1

        xwt[1].append(lat_wt)

        lon_ix1 = np.where(lon < track[i,2])[0]
        if lon_ix1.shape[0] == 0:
            lon_ix1 = 0
        else:
            lon_ix1 = lon_ix1[-1]
        lon_ix2 = lon_ix1 + 1

        if lon_ix2 == lon.shape[0]:
            if verbose:
                sys.stdout.write('NOTE: longitude is unbounded, giving all averaging weight to nearest observation\n')
            lon_ix1 -= 1
            lon_ix2 -= 1
            lon_wt = 0
        else:
            dx1 = lon[lon_ix2] - lon[lon_ix1]
            dx2 = track[i,2] - lon[lon_ix1]
            lon_wt = (dx1 - dx2) / dx1

        xwt[2].append(lon_wt)

        # ---------------------------------------------------------------------
        # Interpolation part 1 - get bounding profiles, interp over time
        # ---------------------------------------------------------------------
        D3 = wt*data[t1,:,lat_ix1:lat_ix2+1,lon_ix1:lon_ix2+1] + (1-wt)*data[t2,:,lat_ix1:lat_ix2+1,lon_ix1:lon_ix2+1]

        # ---------------------------------------------------------------------
        # Interpolation part 2 - interp over lat + lon
        # ---------------------------------------------------------------------
        if np.isnan(D3.flatten()).any():
            if verbose:
                sys.stdout.write('Bounding climatological profile(s) missing data')
                sys.stdout.write(' - taking simple average of available data.\n')
            woa_interp[:,i] = np.nanmean(np.nanmean(D3, axis=2), axis=1)
        else:
            D2 = lat_wt*D3[:,0,:] + (1 - lat_wt)*D3[:,1,:]
            woa_interp[:,i] = lon_wt*D2[:,0] + (1 - lon_wt)*D2[:,1]

    xwt = np.array(xwt)

    return woa_interp, xwt, yrday