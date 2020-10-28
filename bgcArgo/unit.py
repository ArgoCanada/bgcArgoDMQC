import warnings

import numpy as np

# soft attempt to load gsw, but allow for seawater as well
try: 
    import gsw
    flagSA = True
except:
    try:
        # if this also fails, just load gsw to throw the error
        from seawater import pden
        flagSA = False
        warnings.warn('gsw package for thermodynamic equations of seawater not installed, attempting to load seawater package, however seawater is deprecated in favour of gsw-python, see https://teos-10.github.io/GSW-Python/\n')
    except:
        import gsw

def oxy_sol(S, T, unit='micromole/kg'):
    '''
    Calculate oxygen saturation concentration in seawater as a function of
    S & T, in equilibrium with standard coponsition moist air at 1atm total
    pressure. From Garcia & Gordon (1992) eq. 8 (p. 1310) using coefficients
    of Benson & Krause in table 1, as used in Sarmiento & Gruber's "Ocean
    Biogeochemical Dynamics" ch. 3, p. 81, table 3.2.4.
    
    INPUT:
              S: salinity, psu
              T: temperature, deg C
              unit: micromole/kg or millimole/m3, default if micromole/kg
    
    OUTPUT:
              O2sol: oxygen solubility, unit same as input unit
    
    AUTHOR:   Christopher Gordon
              Fisheries and Oceans Canada
              chris.gordon@dfo-mpo.gc.ca
    
    ACKNOWLEDGEMENT: this code is adapted from the SOCCOM SAGE_O2Argo matlab
    code, available via https://github.com/SOCCOM-BGCArgo/ARGO_PROCESSING,
    written by Josh Plant
    
    LAST UPDATE: 21-04-2020
    
    CHANGE LOG:
    '''

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

def pH2O(T, S=0, unit='Pa'):
    '''
    Calculate vapor pressure of water
    
    INPUT:
              T: temperature, deg C
              S: salinity, only necessary for mbar unit
    
    OUTPUT:
              vapor_pressure: vapor pressure of water, Pa or mbar
    
    AUTHOR:   Christopher Gordon
              Fisheries and Oceans Canada
              chris.gordon@dfo-mpo.gc.ca
    
    LAST UPDATE: 06-05-2020
    
    CHANGE LOG:
    '''
    
    # temperature in kelvin
    Tk = T + 273.15

    if unit == 'Pa':
        # from Johnson et al. (2015)
        vapor_pressure = np.exp(52.57 - (6690.9/Tk) - 4.681*np.log(Tk))
    elif unit == 'mbar':
        # SCOR WG 142
        vapor_pressure = 1013.25 * (np.exp(24.4543 - (67.4509*(100/Tk))) - 4.8489*np.log(((Tk/100)) - 0.000544*S))

    return vapor_pressure

def pO2(DOXY, S, T):

    # temperature in kelvin
    Tk = T + 273.15
    # scaled temperature
    Ts = np.log((298.15 - T)/Tk)

    # Benson & Krause (cm^3/dm^3)
    # temperature coefficients
    pA = [3.88767, -0.256847, 4.94457, 4.05010, 3.22014, 2.00907]
    # salinity coefficients
    pB = [-8.17083e-3, -1.03410e-2, -7.37614e-3, -6.24523e-3]
    Co = -4.88682e-7
    O2_volume =  22.3916
    XO2 = 0.20946

    # Used for conversions from [O2] to pO2:
    p1  = np.polyval(pB,Ts)
    S_corr = S*p1 + Co*S**2
    L      = np.polyval(pA,Ts) + S_corr;
    # Oxygen solubility real gas, mmol/m3
    O2sol  = (1000/O2_volume) * np.exp(L)

    PPOX_DOXY = DOXY/O2sol * (1013.25 - pH2O(T)) * XO2

    return PPOX_DOXY

def atmos_pO2(P, pH2O):
    # molar fraction of oxygen in air
    XO2 = 0.20946
    # reference partial pressure of oxygen in air
    pO2 = (P - pH2O) * XO2

    return pO2

# -----------------------------------------------------------------------------
# Section - conversion code from "SCOR WG 142: Quality Control Procedures for 
# Oxygen and Other Biogeochemical Sensors on Floats and Gliders". Code
# adapted from matlab code by Henry Bittig. https://doi.org/10.13155/45915
# -----------------------------------------------------------------------------

def doxy_to_pO2(O2conc, S, T, P=0):
    '''
    convert molar oxygen concentration to oxygen partial pressure

    inputs:
        O2conc - oxygen concentration in umol L-1
        T      - temperature in deg C
        S      - salinity (PSS-78)
        P      - hydrostatic pressure in dbar (default: 0 dbar)

    output:
        pO2    - oxygen partial pressure in mbar

    according to recommendations by SCOR WG 142 "Quality Control Procedures
    for Oxygen and Other Biogeochemical Sensors on Floats and Gliders"

    Written in matlab by: 
    Henry Bittig
    Laboratoire d'Oceanographie de Villefranche-sur-Mer, France
    bittig@obs-vlfr.fr
    28.10.2015
    19.04.2018, v1.1, fixed typo in B2 exponent

    Translated to python by:
    Christopher Gordon
    Bedford Institute of Oceanography, Fisheries and Oceans Canada
    chris.gordon@dfo-mpo.gc.ca
    02.10.2020
    '''

    Tk      = T + 273.15 # temperature in kelvin
    xO2     = 0.20946 # mole fraction of O2 in dry air (Glueckauf 1951)
    pH2Osat = pH2O(T, S=S, unit='mbar') # saturated water vapor in mbar
    # scaled temperature for use in TCorr and SCorr
    sca_T   = np.log((298.15 - T)/(Tk))
    # temperature correction part from Garcia and Gordon (1992), Benson and Krause (1984) refit mL(STP) L-1; and conversion from mL(STP) L-1 to umol L-1
    Tcorr   = 44.6596 * np.exp(2.00907 + 3.22014*sca_T + 4.05010*sca_T**2 + 4.94457*sca_T**3 - 2.56847e-1*sca_T**4 + 3.88767*sca_T**5)
    # salinity correction part from Garcia and Gordon (1992), Benson and Krause (1984) refit ml(STP) L-1
    Scorr   = np.exp(S*(-6.24523e-3 - 7.37614e-3*sca_T - 1.03410e-2*sca_T**2 - 8.17083e-3*sca_T**3) - 4.88682e-7*S**2)
    Vm      = 0.317 # molar volume of O2 in m3 mol-1 Pa dbar-1 (Enns et al. 1965)
    R       = 8.314 # universal gas constant in J mol-1 K-1

    pO2 = O2conc*(xO2*(1013.25 - pH2Osat))/(Tcorr*Scorr)*np.exp(Vm*P/(R*Tk))

    return pO2

def pO2_to_doxy(pO2, S, T, P=0):
    '''
    convert oxygen partial pressure to molar oxygen concentration

    inputs:
        pO2    - oxygen partial pressure in mbar
        T      - temperature in deg C
        S      - salinity (PSS-78)
        P      - hydrostatic pressure in dbar (default: 0 dbar)

    output:
        DOXY   - oxygen concentration in umol L-1

    according to recommendations by SCOR WG 142 "Quality Control Procedures
    for Oxygen and Other Biogeochemical Sensors on Floats and Gliders"


    Written in matlab by: 
    Henry Bittig
    Laboratoire d'Oceanographie de Villefranche-sur-Mer, France
    bittig@obs-vlfr.fr
    28.10.2015
    19.04.2018, v1.1, fixed typo in B2 exponent

    Translated to python by:
    Christopher Gordon
    Bedford Institute of Oceanography, Fisheries and Oceans Canada
    chris.gordon@dfo-mpo.gc.ca
    02.10.2020
    '''

    Tk      = T + 273.15 # temperature in kelvin
    xO2     = 0.20946 # mole fraction of O2 in dry air (Glueckauf 1951)
    pH2Osat = pH2O(T, S=S, unit='mbar') # saturated water vapor in mbar
    # scaled temperature for use in TCorr and SCorr
    sca_T   = np.log((298.15 - T)/(Tk))
    # temperature correction part from Garcia and Gordon (1992), Benson and Krause (1984) refit mL(STP) L-1; and conversion from mL(STP) L-1 to umol L-1
    Tcorr   = 44.6596 * np.exp(2.00907 + 3.22014*sca_T + 4.05010*sca_T**2 + 4.94457*sca_T**3 - 2.56847e-1*sca_T**4 + 3.88767*sca_T**5)
    # salinity correction part from Garcia and Gordon (1992), Benson and Krause (1984) refit ml(STP) L-1
    Scorr   = np.exp(S*(-6.24523e-3 - 7.37614e-3*sca_T - 1.03410e-2*sca_T**2 - 8.17083e-3*sca_T**3) - 4.88682e-7*S**2)
    Vm      = 0.317 # molar volume of O2 in m3 mol-1 Pa dbar-1 (Enns et al. 1965)
    R       = 8.314 # universal gas constant in J mol-1 K-1

    O2conc  = pO2/(xO2*(1013.25 - pH2Osat))/(Tcorr*Scorr)*np.exp(Vm*P/(R*Tk))

    return O2conc

def umol_per_sw_to_mmol_per_L(doxy, S, T, P, Pref=0, lat=None, lon=None):
    '''
    Convert dissolved oxygen concentration in umol kg-1 to mmol L-1.

    Args:
        doxy (float or array-like): dissolved oxygen in umol kg-1
        S (float or array-like): salinity, array of same length as `doxy` or single value
        T (float or array-like): temperature (deg C), array of same length  as `doxy` or single value
        P (float or array-like): pressure (dbar), array of same length  as `doxy` or single value
        Pref (optional, float): reference pressure (dbar) for potential density calculation, default 0
        lat (optional, float or array-like): latitude (deg) for absolute salinity calculation, optional but highly encouraged, function will use practical salinity and produce warning without it
        lon (optional, float or array-like): longitude (deg) for absolute salinity calculation, optional but highly encouraged, function will use practical salinity and produce warning without it

    Returns:
        mmol_L_conc (float or array-like): dissolved oxygen concentration in mmol L-1
    '''

    if flagSA:
        if lat is None and lon is None:
            pot_density = gsw.pot_rho_t_exact(S, T, P)
            warnings.warn('No coordinate information required, proceeding with calculation using practical salinity instead of absolute salinity')
        else:
            pot_density = gsw.pot_rho_t_exact(gsw.SA_from_SP(S, P, lon, lat), T, P)
    else:
        pot_density = pden(S, T, P)

    mmol_L_conc = 1000*doxy / pot_density

    return mmol_L_conc