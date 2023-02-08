#!/usr/bin/python

import warnings

from pathlib import Path
from netCDF4 import Dataset

import unittest
import numpy as np
import pandas as pd
import bgcArgoDMQC as bgc

class coreTest(unittest.TestCase):

    def setUp(self):
        warnings.filterwarnings(action='ignore')

        bgc.set_dirs(
            argo_path=Path(__file__).absolute().parent / 'test_data/Argo/dac',
            woa_path=Path(__file__).absolute().parent / 'test_data/WOA18',
            ncep_path=Path(__file__).absolute().parent / 'test_data/NCEP',
        )

    def test_index_files(self):

        bgc.io.check_index(mode='install')
        bgc.io.check_index()
        bgc.io.update_index()
        
        # get index files
        bgc_index  = bgc.get_index('bgc')
        core_index = bgc.get_index('global')
        syn_index  = bgc.get_index('synthetic')

        self.assertIs(type(bgc_index), pd.core.frame.DataFrame)
        self.assertIs(type(core_index), pd.core.frame.DataFrame)
        self.assertIs(type(syn_index), pd.core.frame.DataFrame) 

    def test_response_time_correction(self):
        # time correction
        time = np.arange(0,20,1)
        temp = 12*np.random.rand(time.shape[0])
        doxy = 200*np.random.rand(time.shape[0])

        doxy_adj = bgc.correct_response_time(time, doxy, temp, 120)
        self.assertIs(type(doxy_adj), np.ndarray)

        doxy_adj_Tconst = bgc.correct_response_time_Tconst(time, doxy, 70)
        self.assertIs(type(doxy_adj_Tconst), np.ndarray)

    def test_pO2(self):
        time = np.arange(0,20,1)
        psal = 35*np.random.rand(time.shape[0])
        temp = 12*np.random.rand(time.shape[0])
        doxy = 200*np.random.rand(time.shape[0])
        pO2 = bgc.unit.pO2(doxy, psal, temp)
        self.assertIs(type(pO2), np.ndarray)

    def test_unit_conversion(self):
        S = np.random.rand(20)
        T = np.random.rand(20)
        P = np.random.rand(20)

        sol = bgc.unit.oxy_sol(S, T, P)

        self.assertIs(type(sol), np.ndarray)

    def test_SCOR_conversion(self):
        doxy = 150 + 50*np.random.rand(30)
        S = 35 + 1.5*np.random.rand(30)
        T = 15 + 5*np.random.rand(30)
        pO2 = bgc.unit.doxy_to_pO2(doxy, S, T)
        doxy_back = bgc.unit.pO2_to_doxy(pO2, S, T)

        self.assertIs(type(pO2), np.ndarray)
        self.assertIs(type(doxy_back), np.ndarray)
