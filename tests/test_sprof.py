#!/usr/bin/python

from pathlib import Path

from netCDF4 import Dataset

import matplotlib.pyplot as plt

import unittest
import numpy as np
import pandas as pd
import bgcArgoDMQC as bgc

class sprofTest(unittest.TestCase):

    def setUp(self):

        print(bgc.resource.path('Argo'))

        bgc.set_dirs(
            argo_path=bgc.resource.path('Argo'),
            ncep_path=bgc.resource.path('NCEP'),
            woa_path=bgc.resource.path('WOA18')
        )

    def test_sprof(self):
        wmo = 4901784
        sprof = bgc.sprof(wmo)
        sprof.clean()
        sprof.reset()
        sprof.clean(bad_flags=[3,4])
        df = sprof.to_dataframe()
        sprof.reset()

        sprof.describe()

        self.assertIsInstance(sprof, bgc.sprof)
        self.assertIs(type(df), pd.core.frame.DataFrame)

    def test_calc_gains(self):
        wmo = 4901784
        sprof = bgc.sprof(wmo)
        ncep_gains = sprof.calc_gains()
        # woa_gains  = sprof.calc_gains(ref='WOA')

        self.assertIs(type(ncep_gains), np.ndarray)
        # self.assertIs(type(woa_gains), np.ndarray)

if __name__ == '__main__':

    unittest.main()