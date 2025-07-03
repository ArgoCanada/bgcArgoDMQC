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
        cyc = 208

        # information for HISTORY_ variables
        history = {
            "HISTORY_INSTITUTION":"BI",
            "HISTORY_STEP":"ARGQ",
            "HISTORY_ACTION":"CF",
            "HISTORY_PARAMETER":"DOXY",
        }

        # load profile
        prof = bgc.prof(wmo=wmo, cycle=cyc, kind='B')
        # load again so that file option can be used
        prof = bgc.prof(file=f'test_data/Argo/dac/meds/{wmo}/profiles/BR{wmo}_{cyc:02d}.nc')

        # clean and reset
        prof.clean(bad_flags=[3,4])
        prof.reset()

        # update DOXY_QC values
        print(prof.DOXY_QC)
        prof.update_field('DOXY_QC', 3, where=prof.DOXY_QC == 1)
        print(prof.DOXY_QC)
        # hard reset
        prof.reset(hard=True)

        # extract dataframe, dict
        df = prof.to_dataframe()
        prof_dict = prof.to_dict()

        self.assertIsInstance(prof, bgc.prof)
        self.assertIs(type(df), pd.core.frame.DataFrame)
        self.assertIs(type(prof_dict), dict)
