#!/usr/bin/python

import numpy as np
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

        # define a function to create an oxygen profile with arbitrary properties
        def oxygen_profile(depth, oxygen_range, central_value, max_grad, depth_grad):

            oxygen = oxygen_range/2*np.tanh((-2*max_grad/oxygen_range)*(depth - depth_grad)) + central_value

            return oxygen

        # just as an example, show the truth, sampled, and corrected
        depth = np.arange(0, 2000, 0.5)[::-1]
        oxy_range = 200 # umol kg-1
        max_grad = 7 # umol kg-1 dbar-1
        oxycline = 100 # dbar

        wmo = 4901784
        syn = bgc.sprof(wmo)

        # create a fake DOXY profile, call it independent
        doxy_mean = syn.DOXY.mean()
        oxygen = oxygen_profile(depth, oxy_range, doxy_mean, max_grad, oxycline)

        data = dict(
            PRES=depth,
            DOXY=oxygen
        )

        syn.add_independent_data(syn.SDN[1], lat=syn.LATITUDE[1], lon=syn.LONGITUDE[1], label='Fake Oxygen', data_dict=data)
        
        with self.assertRaises(ValueError):
            syn.add_independent_data(syn.SDN[1], lat=syn.LATITUDE[1], lon=syn.LONGITUDE[1], label='Fake Oxygen', data_dict=data, TEMP=2)
        
        g = syn.compare_independent_data()
        self.assertIsInstance(g, bgc.plot.pltClass)
