#!/usr/bin/python

from pathlib import Path

from netCDF4 import Dataset

import unittest
import numpy as np
import pandas as pd
import bgcArgoDMQC as bgc

class coreTest(unittest.TestCase):

    tmp = Path('./tmp')
    bgc.set_dirs(
        argo_path=tmp / 'Argo/dac',
        ncep_path=tmp / 'NCEP',
        woa_path=tmp / 'WOA18'
    )

    bgc.io.check_index(mode='install')
    bgc.io.check_index()
    bgc.io.update_index()

    def test_index_files(self):
        # get index files
        bgc_index  = bgc.get_index('bgc')
        core_index = bgc.get_index('global')
        syn_index  = bgc.get_index('synthetic')

        self.assertIs(type(bgc_index), pd.core.frame.DataFrame)
        self.assertIs(type(core_index), pd.core.frame.DataFrame)
        self.assertIs(type(syn_index), pd.core.frame.DataFrame)

    def test_qc_read(self):
        # read QC test
        nc = Dataset(Path('tmp/Argo/dac/meds/4901784/profiles/BD4901784_001.nc'))
        qcp, qcf = bgc.read_history_qctest(nc)
        bgc.util.display_qctests(qcp, qcf)

        self.assertIs(type(qcp), str)
        self.assertIs(type(qcf), str)

    def test_information_criteria(self):
        # aic and bic
        data  = np.random.randn(30)
        resid = np.random.randn(10)/10
        aic = bgc.util.aic(data, resid)

        self.assertIs(type(aic), float)

        data  = np.random.randn(30)
        resid = np.random.randn(10)/10
        bic = bgc.util.bic(data, resid)

        self.assertIs(type(bic), float)   

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
        syn = bgc.sprof(4901784)
        pO2 = bgc.unit.pO2(syn.DOXY, syn.PSAL, syn.TEMP)
        self.assertIs(type(pO2), np.ndarray)

    def test_read_gain_value(self):
        g, eq, comment = bgc.util.read_gain_value(Path('tmp/Argo/dac/meds/4901784/profiles/BD4901784_001.nc'))

        self.assertIs(type(g[0]), np.str_)
        self.assertIs(type(eq[0]), np.str_)
        self.assertIs(type(comment[0]), np.str_)

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

if __name__ == '__main__':
    unittest.main()
