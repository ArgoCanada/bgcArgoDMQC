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
        bgc_index  = bgc.get_index('bgc-b', wmo=4901784)
        core_index = bgc.get_index('core')
        syn_index  = bgc.get_index('bgc-s')
        meta_index  = bgc.get_index('meta', source='remote')
        traj_index  = bgc.get_index('traj', source='local')

        with self.assertRaises(ValueError):
            bgc.get_index('not a real index name')

        self.assertIs(type(bgc_index), pd.core.frame.DataFrame)
        self.assertIs(type(core_index), pd.core.frame.DataFrame)
        self.assertIs(type(syn_index), pd.core.frame.DataFrame) 
        self.assertIs(type(meta_index), pd.core.frame.DataFrame) 
        self.assertIs(type(traj_index), pd.core.frame.DataFrame) 

    def test_response_time_correction(self):
        # time correction
        time = np.arange(0, 5000, 10)
        temp = 12*np.random.rand(time.shape[0])
        def oxygen_profile(depth, oxygen_range, central_value, max_grad, depth_grad):

            oxygen = oxygen_range/2*np.tanh((-2*max_grad/oxygen_range)*(depth - depth_grad)) + central_value

            return oxygen

        # just as an example, show the truth, sampled, and corrected
        depth = np.arange(0, 500, 1)[::-1]
        oxy_range = 200 # umol kg-1
        central = 120 # umol kg-1
        max_grad = 7 # umol kg-1 dbar-1
        oxycline = 100 # dbar
        oxygen = oxygen_profile(depth, oxy_range, central, max_grad, oxycline)

        vel = np.mean(np.diff(depth)/np.diff(time))
        Il = bgc.estimate_boundary_layer(vel)

        doxy_adj = bgc.correct_response_time(time, oxygen, temp, Il)
        self.assertIs(type(doxy_adj), np.ndarray)

        doxy_adj_Tconst = bgc.correct_response_time_Tconst(time, oxygen, 70)
        self.assertIs(type(doxy_adj_Tconst), np.ndarray)

        doxy_original = bgc.sample(time, doxy_adj_Tconst, temp, Il)
        self.assertTrue(np.all(np.abs(doxy_original - oxygen) < 1))

        doxy_original_Tconst = bgc.sample_Tconst(time, doxy_adj_Tconst, 70)
        self.assertTrue(np.all(np.abs(doxy_original_Tconst - oxygen) < 1))

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
        ml_per_l = bgc.unit.mL_per_L_to_umol_per_L(np.array(10*[6]), np.array(10*[13]))
        pH2O_mbar = bgc.unit.pH2O(T, unit='mbar')
        pH2O_Pa = bgc.unit.pH2O(T, unit='mbar')
        
        with self.assertRaises(ValueError):
            pH2O_Pa = bgc.unit.pH2O(T, unit='not a real unit')

        self.assertIs(type(pO2), np.ndarray)
        self.assertIs(type(doxy_back), np.ndarray)
        self.assertIs(type(ml_per_l), np.ndarray)
        self.assertIs(type(pH2O_mbar), np.ndarray)
        self.assertIs(type(pH2O_Pa), np.ndarray)
    
    def gridded_var_read(self):
        
        wmo = 6902870
        argo_path = Path(__file__).absolute().parent / 'test_data/Argo/dac'
        argo_path.mkdir(exist_ok=True)
        bgc.io.get_argo(wmo, local_path=argo_path, overwrite=False, nfiles=5)

        nc = Dataset(argo_path / f'coriolis/{wmo}/{wmo}_Sprof.nc')
        float_dict = bgc.read_gridded_variables(nc)

        self.assertIs(type(float_dict), dict)
        self.assertIs(type(float_dict['TEMP']), np.ndarray)

