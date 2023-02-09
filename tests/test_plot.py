#!/usr/bin/python

import matplotlib.pyplot as plt
import unittest
from pathlib import Path

import bgcArgoDMQC as bgc

class plottingTest(unittest.TestCase):

    def setUp(self):

        bgc.set_dirs(
            argo_path=Path(__file__).absolute().parent / 'test_data/Argo/dac',
            woa_path=Path(__file__).absolute().parent / 'test_data/WOA18',
            ncep_path=Path(__file__).absolute().parent / 'test_data/NCEP',
        )

        wmo = 4901784
        syn = bgc.sprof(wmo)

    def test_profile_plot(self):
        wmo = 4901784
        syn = bgc.sprof(wmo)

        g_pres = syn.plot('qcprofiles', varlist=['TEMP', 'DOXY'])
        self.assertIsInstance(g_pres, bgc.plot.pltClass)
        plt.close(g_pres.fig)

        g_pden = syn.plot('profiles', varlist=['PSAL', 'DOXY'], Nprof=5, Ncycle=3, zvar='PDEN')
        self.assertIsInstance(g_pden, bgc.plot.pltClass)
        plt.close(g_pden.fig)

    def test_scatter_plot(self):

        wmo = 4901784
        syn = bgc.sprof(wmo)

        g = syn.plot('cscatter', varname='DOXY')
        self.assertIsInstance(g, bgc.plot.pltClass)
        plt.close(g.fig)

    def test_independent_data(self):

        return
