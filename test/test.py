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

if __name__ == '__main__':
    unittest.main()