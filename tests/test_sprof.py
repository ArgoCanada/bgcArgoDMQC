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

        bgc.set_dirs(
            argo_path=Path('test_data/Argo/dac').absolute(),
            ncep_path=Path('test_data/NCEP').absolute(),
            woa_path=Path('test_data/WOA18').absolute()
        )

        wmo = 4901784
        bgc.io.Path.ARGO_PATH.mkdir(exist_ok=True, parents=True)
        bgc.io.get_argo(wmo, local_path=bgc.io.Path.ARGO_PATH, overwrite=False, nfiles=2)

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
