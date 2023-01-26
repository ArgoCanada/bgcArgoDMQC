import numpy as np
from scipy.interpolate import interp1d, interp2d

def oxy_b(dt, tau):
    '''
    Calculates the coefficient _b_ from equation 4 of Bittig et al. 2014

    Args:
        - dt: change in time (s)
        - tau: time constant (s)
    Returns: 
        - 1/(1 + 2(tau/dt))
    '''

    inv_b = 1 + 2*(tau/dt)
    return 1/inv_b

def oxy_a(dt, tau):
    '''
    Calculates the coefficient _a_ from equation 4 of Bittig et al. 2014
    
    Args:
        - dt: change in time (s)
        - tau: time constant (s)
    Returns: 
        - 1 - 2b = 1 - 2/(1 + 2(tau/dt))
    '''
    return 1 - 2*oxy_b(dt, tau)

def correct_response_time_Tconst(t, O2, tau):
    '''
    Corrects optode time response error while assuming constant temperature

    Args:
        - t: time in days (does not have to begin with 0, could take input directly from an Argo netCDF 'JULD' variable for example)
        - O2: oxygen, can be in concentration or partial pressure units
        - tau: time constant in seconds
    Returns:
        - corrected oxygen on the same time axis
    '''

    # array for the loop
    N = O2.shape[0]
    mean_oxy  = np.array((N-1)*[np.nan])
    mean_time = np.array((N-1)*[np.nan])

    # convert time to seconds
    t_sec = t*24*60*60

    # loop through oxygen data
    for i in range(N-1):
        dt = t_sec[i+1] - t_sec[i]

        # do the correction using the mean filter, get the mean time
        mean_oxy[i]  = (1/(2*oxy_b(dt, tau)))*(O2[i+1] - oxy_a(dt, tau)*O2[i])
        mean_time[i] = t_sec[i] + dt/2
    
    # interpolate back to original times for output
    f = interp1d(mean_time, mean_oxy, kind='linear', bounds_error=False, fill_value='extrapolate')
    O2_out = f(t_sec)

    return O2_out

from ..lut import lut as lut_data

def correct_response_time(t, O2, T, thickness):
    '''
    Corrects optode response time error using a temperature dependent time
    constant calculated using the lookup table provided in Bittig and
    Kortzinger 2017.

    Args:
        - t: time in days (does not have to begin with 0, could take input directly from an Argo netCDF 'JULD' variable for example)
        - O2: oxygen, can be in concentration or partial pressure units
        - T: temperature, deg C
        - thickness: boundary layer thickness, micrometers
    Returns:
        - corrected oxygen on the same time axis
    '''

    # convert time to seconds
    t_sec = t*24*60*60

    # array for the loop
    N = O2.shape[0]
    mean_oxy  = np.array((N-1)*[np.nan])
    mean_time = t_sec[:-1] + np.diff(t_sec)/2
    mean_temp = T[:-1] + np.diff(T)/2

    tau_T = lookup_tau(thickness, mean_temp)

    # loop through oxygen data 
    for i in range(N-1):
        dt = t_sec[i+1] - t_sec[i]

        # do the correction using the mean filter, get the mean time
        mean_oxy[i]  = (1/(2*oxy_b(dt, tau_T[i])))*(O2[i+1] - oxy_a(dt, tau_T[i])*O2[i])
    
    # interpolate back to original times for output
    f = interp1d(mean_time, mean_oxy, kind='linear', bounds_error=False, fill_value='extrapolate')
    O2_out = f(t_sec)

    return O2_out

def sample_Tconst(t, O2, tau):
    '''
    Simulate how an oxygen optode samples an oxygen profile with time response
    error using a low pass filter, assuming constant temperature

    Args:
        - t: time in days (does not have to begin with 0, could take input directly from an Argo netCDF 'JULD' variable for example)
        - O2: oxygen profile to be sampled, can be in concentration or partial pressure units
        - tau: time constant in seconds
    Returns:
        - sampled/filtered oxygen, at the same vertical resolution as the input oxygen
    '''
    # convert time to seconds
    t_sec = t*24*60*60

    N = O2.shape[0]
    c_filt = np.nan*np.ones((N,))
    c_filt[0] = O2[0]

    for i in range(N-1):
        dt = t_sec[i+1] - t_sec[i]
        c_filt[i+1] = oxy_a(dt, tau)*c_filt[i] + oxy_b(dt, tau)*(O2[i+1]+O2[i])
    
    return c_filt

def sample(t, O2, T, thickness):
    '''
    Simulate how an oxygen optode samples an oxygen profile with time response 
    error using a low pass filter, using a temperature dependent time constant 
    for a given boundary layer thickness

    Args:
        - t: time in days (does not have to begin with 0, could take input directly from an Argo netCDF 'JULD' variable for example)
        - O2: oxygen profile to be sampled, can be in concentration or partial pressure units
        - T: temperature, deg C
        - thickness: boundary layer thickness, micrometers
    Returns:
        - sampled/filtered oxygen, at the same vertical resolution as the input oxygen
    '''

    # convert time to seconds
    t_sec = t*24*60*60

    mean_temp = T[:-1] + np.diff(T)/2

    N = O2.shape[0]
    c_filt = np.nan*np.ones((N,))
    c_filt[0] = O2[0]

    tau_T = lookup_tau(thickness, mean_temp)

    # loop through oxygen data 
    for i in range(N-1):
        dt = t_sec[i+1] - t_sec[i]

        c_filt[i+1] = oxy_a(dt, tau_T[i])*c_filt[i] + oxy_b(dt, tau_T[i])*(O2[i+1]+O2[i])
    
    return c_filt

def lookup_tau(thickness, T):
    '''
    Function to lookup the time constant (tau) value for a given boundary layer 
    thickness and temperature

    Args: 
        - T: temperature, deg C
        - thickness: boundary layer thickness, micrometers
    Returns: 
        - time constant, seconds    
    '''
    # load temperature, boundary layer thickness, and tau matrix from 
    # look-up table provided in the supplement to Bittig and Kortzinger (2017)
    N = T.shape[0]
    lut_lL = lut_data[0,1:]
    lut_T  = lut_data[1:,0]
    tau100 = lut_data[1:,1:]
    thickness = thickness*np.ones((N,))

    # translate boundary layer thickness to temperature dependent tau
    f_thickness = interp2d(lut_T, lut_lL, tau100.T, bounds_error=False)
    tau_T = np.squeeze(f_thickness(T, thickness))[0,:]

    return tau_T

def estimate_boundary_layer(vel):
    """
    Estimate boundary layer thickness based on platform velocity. Follows eq.
    17 from Bittig et al. (2018).
    """

    if vel <= 0.095:
        return 210 - (110/0.095)*np.abs(vel)
    else:
        return 20 + (80/0.905)*(1 - np.abs(vel))
