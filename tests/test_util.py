#!/usr/bin/python

import warnings
from pathlib import Path

import unittest

from netCDF4 import Dataset
import numpy as np
import bgcArgoDMQC as bgc

class coreTest(unittest.TestCase):

    def setUp(self):
        warnings.filterwarnings(action='ignore')

        bgc.set_dirs(
            argo_path=Path(__file__).absolute().parent / 'test_data/Argo/dac',
            woa_path=Path(__file__).absolute().parent / 'test_data/WOA18',
            ncep_path=Path(__file__).absolute().parent / 'test_data/NCEP',
        )
    
    def test_qc_read(self):
        # read QC test
        nc = Dataset(bgc.io.Path.ARGO_PATH / 'meds/4901784/profiles/BD4901784_001.nc')
        qcp, qcf = bgc.read_history_qctest(nc)
        bgc.util.display_qctests(qcp, qcf)

        self.assertIs(type(qcp), np.str_)
        self.assertIs(type(qcf), np.str_)
    
    def test_information_criteria(self):
        # aic and bic
        data  = np.random.randn(30)
        resid = np.random.randn(150)
        aic = bgc.util.aic(data, resid)

        self.assertIs(type(aic), np.float64)

        data  = np.random.randn(30)
        resid = np.random.randn(150)
        bic = bgc.util.bic(data, resid)

        self.assertIs(type(bic), np.float64) 

    def test_haversine(self):

        c1 = (45, -60)
        c2 = (46, -61)

        distance = bgc.util.haversine(c1, c2)
        self.assertGreaterEqual(distance, 0)

    def test_file_vars(self):
        files = list((bgc.io.Path.ARGO_PATH / 'meds/4901784/profiles/').glob('*.nc'))
        varnames = bgc.util.get_vars(files)

        self.assertIs(type(varnames), list)

    def test_worst_flag(self):

        #     0  1  2  3  4  5
        f1 = [1, 3, 1, 2, 4, 4]
        f2 = [2, 9, 1, 5, 3, 3]
        #     2  9  1  5  4  4

        ff = bgc.util.get_worst_flag(np.array(f1), np.array(f2))
        self.assertEqual(ff[0], 2)
        self.assertEqual(ff[5], 4)
        self.assertEqual(ff[1], 9)

    def test_nc_array_functions(self):

        N = 2
        M = N+1

        for i in range(4):
            dim = i+1
            for ax in range(dim):
                if dim == 1:
                    old_array = np.array(N*[1])
                    new_array = np.array(M*[2])
                    final = bgc.util.refill_array(ax, dim, new_array, old_array)
                    self.assertTrue(np.all(final[:N] == 1))
                    self.assertTrue(final[-1] == 2)
                elif dim == 2:
                    old_array = np.ones((N, N))
                    if ax == 0:
                        new_array = 2*np.ones((M, N))
                    else:
                        new_array = 2*np.ones((N, M))
                    final = bgc.util.refill_array(ax, dim, new_array, old_array)
                    final = bgc.util.refill_array(ax, dim, new_array, old_array)
                    self.assertTrue(np.all(final[:N, :N] == 1))
                    self.assertTrue(final[-1, -1] == 2)
                elif dim == 3:
                    old_array = np.ones((N, N, N))
                    if ax == 0:
                        new_array = 2*np.ones((M, N, N))
                    elif ax == 1:
                        new_array = 2*np.ones((N, M, N))
                    else:
                        new_array = 2*np.ones((N, N, M))
                    final = bgc.util.refill_array(ax, dim, new_array, old_array)
                    self.assertTrue(np.all(final[:N, :N, :N] == 1))
                    self.assertTrue(final[-1, -1, -1] == 2)
                elif dim == 4:
                    old_array = np.ones((N, N, N, N))
                    if ax == 0:
                        new_array = 2*np.ones((M, N, N, N))
                    elif ax == 1:
                        new_array = 2*np.ones((N, M, N, N))
                    elif ax == 2:
                        new_array = 2*np.ones((N, N, M, N))
                    else:
                        new_array = 2*np.ones((N, N, N, M))
                    final = bgc.util.refill_array(ax, dim, new_array, old_array)
                    self.assertTrue(np.all(final[:N, :N, :N, :N] == 1))
                    self.assertTrue(final[-1, -1, -1, -1] == 2)

    def test_read_gain(self):

        fn = bgc.io.Path.ARGO_PATH / 'meds/4901784/profiles/BD4901784_001.nc'
        bgc.util.read_gain_value(fn)
