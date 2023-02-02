#!/usr/bin/python

import matplotlib.pyplot as plt

import unittest
import bgcArgoDMQC as bgc

class plottingTest(unittest.TestCase):

    def setUp(self):

        wmo = 4901784
        bgc.resource.path('Argo').mkdir(exist_ok=True, parents=True)
        bgc.io.get_argo(wmo, local_path=bgc.resource.path('Argo'), overwrite=False, nfiles=2)
        
        bgc.set_dirs(
            argo_path=bgc.resource.path('Argo'),
            ncep_path=bgc.resource.path('NCEP'),
            woa_path=bgc.resource.path('WOA18')
        )

    def test_gain_plot(self):
        wmo = 4901784
        syn  = bgc.sprof(wmo)

        syn.calc_gains()
        syn.calc_gains(ref='WOA')

        g_ncep = syn.plot('gain', ref='NCEP')
        self.assertIsInstance(g_ncep, bgc.fplt.pltClass)
        plt.close(g_ncep.fig)

        g_woa  = syn.plot('gain', ref='WOA')
        self.assertIsInstance(g_woa, bgc.fplt.pltClass)
        plt.close(g_woa.fig)

    def test_scatter_plot(self):
        wmo = 4901784
        syn = bgc.sprof(wmo)
        g = syn.plot('cscatter', varname='DOXY', ylim=(0,500))

        self.assertIsInstance(g, bgc.fplt.pltClass)

        plt.close(g.fig)

    def test_profile_plot(self):
        wmo = 4901784
        syn = bgc.sprof(wmo)

        g_pres = syn.plot('profiles', varlist=['TEMP', 'DOXY'])
        self.assertIsInstance(g_pres, bgc.fplt.pltClass)
        plt.close(g_pres.fig)

        g_pden = syn.plot('profiles', varlist=['PSAL', 'DOXY'], Nprof=5, Ncycle=3, zvar='PDEN')
        self.assertIsInstance(g_pden, bgc.fplt.pltClass)
        plt.close(g_pden.fig)
