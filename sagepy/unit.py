#!/usr/bin/python

import numpy as np

def oxy_sol(S, T, unit='micromole/kg'):
    # -------------------------------------------------------------------------
    # oxy_sol
    # -------------------------------------------------------------------------
    #
    # Calculate oxygen saturation concentration in seawater as a function of
    # S & T, in equilibrium with standard coponsition moist air at 1atm total
    # pressure. From Garcia & Gordon (1992) eq. 8 (p. 1310) using coefficients
    # of Benson & Krause in table 1, as used in Sarmiento & Gruber's "Ocean
    # Biogeochemical Dynamics" ch. 3, p. 81, table 3.2.4.
    #
    # INPUT:
    #           S: salinity, psu
    #           T: temperature, deg C
    #           unit: micromole/kg or millimole/m3, default if micromole/kg
    #
    # OUTPUT:
    #           O2sol: oxygen solubility, unit same as input unit
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
        raise ValueError('Unrecognized unit string - valid units are ''micromole/kg'' or ''millimole/m3''.')

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

def pH2O(S, T):
    # -------------------------------------------------------------------------
    # pH2O
    # -------------------------------------------------------------------------
    #
    # Calculate vapor pressure of water
    #
    # INPUT:
    #           S: salinity, psu
    #           T: temperature, deg C
    #
    # OUTPUT:
    #           pH2O: vapor pressure of water, Pa
    #
    # AUTHOR:   Christopher Gordon
    #           Fisheries and Oceans Canada
    #           chris.gordon@dfo-mpo.gc.ca
    #
    # LAST UPDATE: 06-05-2020
    #
    # CHANGE LOG:
    #
    # -------------------------------------------------------------------------
    
    # temperature in kelvin
    Tk = T + 273.15

    # define coefficient array
    D = np.array([24.2543, -67.4509, -4.8489, -5.44e-4])
    # compute exponent
    Dsum = D[0] + D[1]*(100/(Tk)) + D[2]*np.log(Tk/100) + D[3]*S

    return 1013.25*np.exp(Dsum)

def pO2(Pncep, pH2O):
    # -------------------------------------------------------------------------
    # pO2
    # -------------------------------------------------------------------------
    #
    # Calculate vapor pressure of water
    #
    # INPUT:
    #           Pncep: NCEP air pressure at sea surface, Pa
    #           pH2O: vapor pressure of water, Pa
    #
    # OUTPUT:
    #           pO2: partial pressure of atmospheric oxygen, Pa
    #
    # AUTHOR:   Christopher Gordon
    #           Fisheries and Oceans Canada
    #           chris.gordon@dfo-mpo.gc.ca
    #
    # LAST UPDATE: 06-05-2020
    #
    # CHANGE LOG:
    #
    # -------------------------------------------------------------------------

    # mole fraction of oxygen 
    XO2 = 0.20946

    return (Pncep - pH2O) * XO2