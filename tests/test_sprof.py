#!/usr/bin/python

from pathlib import Path
import shutil

from netCDF4 import Dataset

import matplotlib.pyplot as plt

import unittest
import numpy as np
import pandas as pd
import bgcArgoDMQC as bgc

global wmo
wmo = 4902481

global data_path
data_path = Path('/Users/gordonc/Documents/data')

bgc.set_dirs(
    argo_path=data_path / 'Argo',
    ncep_path=data_path / 'NCEP',
    woa_path=data_path / 'WOA18'
)

class sprofTest(unittest.TestCase):

    bgc.set_dirs(
        argo_path=data_path / 'Argo',
        ncep_path=data_path / 'NCEP',
        woa_path=data_path / 'WOA18'
    )

    def test_sprof(self):
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
        sprof = bgc.sprof(wmo)
        ncep_gains = sprof.calc_gains()
        woa_gains  = sprof.calc_gains(ref='WOA')

        self.assertIs(type(ncep_gains), np.ndarray)
        self.assertIs(type(woa_gains), np.ndarray)