#!/usr/bin/python

import sys

from pathlib import Path
import ftplib
import gzip

from . import util

def get_woa18(varname, local_path='./', ftype='netcdf'):
    # -------------------------------------------------------------------------
    # get_woa18
    # -------------------------------------------------------------------------
    #
    # Function to load in all data from a single float, using BRtraj, meta,
    # and Sprof files
    #
    # INPUT:
    #           varname: woa18 variable to download data for - one of:
    #               T: temperature
    #               S: salinity
    #               O2: dissolved oxygen
    #               O2sat: oxygen percent saturation
    #               NO3: nitrate
    #               Si: silicate
    #               PO4: phosphate
    #           local_path: path to save files to (inside variable folder, ex.
    #               if local_path = '../data/woa', and varname = 'T', files
    #               will be saved in ../data/woa/temperature/*.nc. Defaults
    #               to current directory
    #           ftype: file format to download, defaults to netcdf (.nc):
    #               ascii: .dat
    #               csv: .csv
    #               netcdf: .nc
    #
    # OUTPUT:
    #
    # AUTHOR:   Christopher Gordon
    #           Fisheries and Oceans Canada
    #           chris.gordon@dfo-mpo.gc.ca
    #
    # LAST UPDATE: 29-04-2020
    #
    # CHANGE LOG:
    #
    # -------------------------------------------------------------------------

    url = 'ftp.nodc.noaa.gov'
    param, ftype, ftpdir = util.decode_woa_var(varname)

    ftp = ftplib.FTP(url, 'anonymous', 'chris.gordon@dfo-mpo.gc.ca')

    

    return ftp

def get_argo(*args):
    # -------------------------------------------------------------------------
    # get_argo
    # -------------------------------------------------------------------------
    #
    # Function to load in all data from a single float, using BRtraj, meta,
    # and Sprof files
    #
    # INPUT:
    #           Inputs may vary depending on desired performance. Multiple
    #           arguments may be provided to download all files from a certain
    #           float or argo defined geographical area. A single path to a 
    #           file may be provided to download that file. A list of files 
    #           may be provided as well.
    #
    # OUTPUT:
    #
    # AUTHOR:   Christopher Gordon
    #           Fisheries and Oceans Canada
    #           chris.gordon@dfo-mpo.gc.ca
    #
    # LAST UPDATE: 29-04-2020
    #
    # CHANGE LOG:
    #
    # -------------------------------------------------------------------------

    return None

# if __name__ == '__main__':
#     sys.stdout.write('Running io.py as __main__ - testing diagnostics of each function\n')

#     ftp = get_woa18('o2sat')