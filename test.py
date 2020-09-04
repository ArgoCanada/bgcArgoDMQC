#!/usr/bin/python

from pathlib import Path

import unittest
import numpy as np
import bgcArgo as bgc

global wmo
wmo = 4902481

bgc.set_dirs(
    argo_path='tmp/Argo',
    ncep_path='tmp/NCEP',
    woa_path='tmp/WOA18'
)

class downloadTest(unittest.TestCase):

    def download_ncep(self):

        self.assertTrue(Path('tmp/NCEP/').exists())

    def download_woa(self):

        self.assertTrue(Path('tmp/WOA18/').exists())

    def download_argo(self):

        for w in [wmo-1, wmo]:

            bgc.io.get_argo(w, local_path='tmp/Argo')

            self.assertTrue(Path('tmp/Argo/{}'.format(w)))

class sprofTest(unittest.TestCase):

    def test_sprof(self):
        sprof = bgc.sprof(wmo)

        self.assertIsInstance(sprof, bgc.sprof)

    def test_calc_gains(self):
        sprof = bgc.sprof(wmo)
        ncep_gains = sprof.calc_gains()
        woa_gains  = sprof.calc_gains(ref='WOA')

        self.assertIs(type(ncep_gains), np.ndarray)
        self.assertIs(type(woa_gains), np.ndarray)

class profilesTest(unittest.TestCase):
    def test_sprof(self):
        prof = bgc.profiles(wmo)
        # test multiple profs
        profs = bgc.profiles([wmo-1, wmo])
        # test specific cycles
        cycs = bgc.profiles(wmo, cycles=np.arange(1,10))

        self.assertIsInstance(prof, bgc.profiles)
        self.assertIsInstance(profs, bgc.profiles)
        self.assertIsInstance(cycs, bgc.profiles)

    def test_calc_gains(self):
        prof = bgc.profiles(wmo)
        woa_gains  = prof.calc_gains(ref='WOA')

        self.assertIs(type(woa_gains), np.ndarray)

if __name__ == '__main__':
    unittest.main()