#!/usr/bin/python

import unittest
import numpy as np
import bgcArgo as bgc

global wmo
wmo = 4902481

bgc.set_dirs(
    argo_path='/Users/gordonc/Documents/data/Argo',
    ncep_path='/Users/gordonc/Documents/data/NCEP',
    woa_path='/Users/gordonc/Documents/data/WOA18'
)

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