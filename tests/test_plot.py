#!/usr/bin/python

import matplotlib.pyplot as plt
import unittest
from pathlib import Path

import bgcArgoDMQC as bgc

class plottingTest(unittest.TestCase):

    def setUp(self):

        bgc.set_dirs(
            argo_path=Path('test_data/Argo/dac').absolute(),
            ncep_path=Path('test_data/NCEP').absolute(),
            woa_path=Path('test_data/WOA18').absolute()
        )

    def test_profile_plot(self):
        wmo = 4901784
        syn = bgc.sprof(wmo)

        g_pres = syn.plot('qcprofiles', varlist=['TEMP', 'DOXY'])
        self.assertIsInstance(g_pres, bgc.plot.pltClass)
        plt.close(g_pres.fig)

        g_pden = syn.plot('profiles', varlist=['PSAL', 'DOXY'], Nprof=5, Ncycle=3, zvar='PDEN')
        self.assertIsInstance(g_pden, bgc.plot.pltClass)
        plt.close(g_pden.fig)

if __name__ == '__main__':
    unittest.main()
