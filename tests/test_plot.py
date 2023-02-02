#!/usr/bin/python

import matplotlib.pyplot as plt

import unittest
import bgcArgoDMQC as bgc

class plottingTest(unittest.TestCase):

    def setUp(self):

        bgc.set_dirs(
            argo_path=bgc.resource.path('Argo'),
            ncep_path=bgc.resource.path('NCEP'),
            woa_path=bgc.resource.path('WOA18')
        )

    def test_profile_plot(self):
        wmo = 4901784
        syn = bgc.sprof(wmo)

        g_pres = syn.plot('profiles', varlist=['TEMP', 'DOXY'])
        self.assertIsInstance(g_pres, bgc.fplt.pltClass)
        plt.close(g_pres.fig)

        g_pden = syn.plot('profiles', varlist=['PSAL', 'DOXY'], Nprof=5, Ncycle=3, zvar='PDEN')
        self.assertIsInstance(g_pden, bgc.fplt.pltClass)
        plt.close(g_pden.fig)
