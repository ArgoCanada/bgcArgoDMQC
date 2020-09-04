#!/usr/bin/python

from pathlib import Path
import shutil

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

    def test_download_ncep(self):

        bgc.io.get_ncep('pres', local_path='tmp/NCEP')
        bgc.io.get_ncep('land', local_path='tmp/NCEP')

        self.assertTrue(Path('tmp/NCEP/pres/pres.sfc.gauss.2010.nc').exists())

    def test_download_woa(self):

        bgc.io.get_woa18('O2sat', local_path='tmp/WOA18')

        self.assertTrue(Path('tmp/WOA18/o2sat/woa18_all_O00_01.nc').exists())

    def test_download_argo(self):

        dac = [bgc.get_dac(w) for w in [wmo-1,wmo]]

        dacpath = '/ifremer/argo/dac'
        fltpath = ['{}/{}/{}'.format(dacpath, d, w) for d, w in zip(dac, [wmo-1,wmo])]
        bgc.io.get_argo(fltpath, local_path='tmp/Argo')

        self.assertTrue(Path('tmp/Argo/{}/{}'.format(dac[0],wmo-1)).exists())
        self.assertTrue(Path('tmp/Argo/{}/{}'.format(dac[1],wmo)).exists())

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
    if not Path('tmp').exists():
        Path('tmp').mkdir()
    if not Path('tmp/Argo').exists():
        Path('tmp/Argo').mkdir()
    if not Path('tmp/NCEP').exists():
        Path('tmp/NCEP').mkdir()
    if not Path('tmp/NCEP/pres').exists():
        Path('tmp/NCEP/pres').mkdir()
    if not Path('tmp/NCEP/land').exists():
        Path('tmp/NCEP/land').mkdir()
    if not Path('tmp/WOA18').exists():
        Path('tmp/WOA18').mkdir()
    if not Path('tmp/WOA18/o2sat').exists():
        Path('tmp/WOA18/o2sat').mkdir()

    unittest.main()
    
    # shutil.rmtree('tmp')