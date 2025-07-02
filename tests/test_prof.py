#!/usr/bin/python

from pathlib import Path

import unittest
import pandas as pd
import bgcArgoDMQC as bgc

class profTest(unittest.TestCase):

    def setUp(self):

        bgc.set_dirs(
            argo_path=Path(__file__).absolute().parent / 'test_data/Argo/dac',
            ncep_path=Path(__file__).absolute().parent / 'test_data/NCEP',
            woa_path=Path(__file__).absolute().parent / 'test_data/WOA18',
        )

        wmo = 4901784
        bgc.io.Path.ARGO_PATH.mkdir(exist_ok=True, parents=True)
        bgc.io.get_argo(wmo, local_path=bgc.io.Path.ARGO_PATH, overwrite=False, nfiles=2)

    def test_prof(self):
        wmo = 4901784
        cyc = 2
        sprof = bgc.prof(wmo, cyc, kind='B')