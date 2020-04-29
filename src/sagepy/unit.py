#!/usr/bin/python

import numpy as np

def oxy_sol(S, T, unit='micromole/kg'):
    # -------------------------------------------------------------------------
    # load_float_data.py
    # -------------------------------------------------------------------------
    #
    # Calculate oxygen saturation concentration in seawater as a function of
    # S & T, in equilibrium with standard coponsition moist air at 1atm total
    # pressure. From Garcia & Gordon (1992) eq. 8 (p. 1310) using coefficients
    # of Benson & Krause in table 1, as used in Sarmiento & Gruber's "Ocean
    # Biogeochemical Dynamics" ch. 3, p. 81, table 3.2.4.
    #
    # INPUT:
    #
    # OUTPUT:
    #
    # AUTHOR:   Christopher Gordon
    #           Fisheries and Oceans Canada
    #           chris.gordon@dfo-mpo.gc.ca
    #
    # ACKNOWLEDGEMENT: this code is adapted from the SOCCOM SAGE_O2Argo matlab
    # code, available via https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING,
    # written by Josh Plant
    #
    # LAST UPDATE: 21-04-2020
    #
    # CHANGE LOG:
    #
    # -------------------------------------------------------------------------

    # check for improper units
    if unit != 'micromole/kg' and unit != 'millimole/m3':
        raise ValueError('Unrecognized unit string - valid units are ''micro'' or ''milli''.')

    if unit == 'micromole/kg':
        A = [3.80369, -9.86643e-2, 5.10006, 4.17887, 3.20291, 5.80871]
        B = [-9.51519e-3, -1.13864e-2, -7.70028e-3, -7.01577e-3]
        C = -2.75915e-7

    elif unit == 'millimole/m3':
        A = [3.88767, -0.256847, 4.94457, 4.05010, 3.22014, 2.00907]
        B = [-8.17083e-3, -1.03410e-2, -7.37614e-3, -6.24523e-3]
        C = -4.88682e-7

    # Scaled temperature
    Ts = np.log((298.15 - T)/(273.15 + T));
    L = np.polyval(A,Ts) + S*np.polyval(B,Ts) + C*S**2

    O2sol = np.exp(L)

    return O2sol
