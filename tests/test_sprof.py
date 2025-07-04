#!/usr/bin/python

from pathlib import Path

import unittest
import pandas as pd
import bgcArgoDMQC as bgc

class sprofTest(unittest.TestCase):

    def setUp(self):

        bgc.set_dirs(
            argo_path=Path(__file__).absolute().parent / 'test_data/Argo/dac',
            ncep_path=Path(__file__).absolute().parent / 'test_data/NCEP',
            woa_path=Path(__file__).absolute().parent / 'test_data/WOA18',
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

        print(sprof.DOXY)
        sprof['DOXY'] = sprof['DOXY']*1.1
        print(sprof['DOXY'])
        sprof.check_range('all')
        sprof_dict = sprof.to_dict()

        sprof.clean(bad_flags=[1, 2])

        sprof.calc_fixed_error()

        with self.assertRaises(ValueError):
            sprof.plot('gain', ref='Not a real ref')
        with self.assertRaises(ValueError):
            sprof.plot('Not a real plot')

        self.assertIsInstance(sprof, bgc.sprof)
        self.assertIs(type(df), pd.core.frame.DataFrame)
