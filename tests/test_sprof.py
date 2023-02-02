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

        wmo = 4901784
        bgc.resource.path('Argo').mkdir(exist_ok=True, parents=True)
        bgc.io.get_argo(wmo, local_path=bgc.resource.path('Argo'), overwrite=False, nfiles=2)

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

if __name__ == '__main__':

    unittest.main()