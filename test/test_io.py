#!/usr/bin/python

from pathlib import Path

import bgcArgoDMQC as bgc

import unittest

class downloadTest(unittest.TestCase):

    def test_download_refdata(self):

        bgc.diagnostic.simple_test(data_path='tmp')
        bgc.io.get_ncep('rhum', local_path='tmp/NCEP', overwrite=True, years=[2019, 2020])

        self.assertTrue(Path('tmp/NCEP/pres/pres.sfc.gauss.2019.nc').exists())
        self.assertTrue(Path('tmp/NCEP/land/land.sfc.gauss.nc').exists())
        self.assertTrue(Path('tmp/NCEP/rhum/rhum.sig995.2019.nc').exists())
        self.assertTrue(Path('tmp/WOA18/o2sat/woa18_all_O00_01.nc').exists())

    def test_download_argo(self):

        bgc.io.get_argo(3900407, local_path='tmp/Argo/aoml', overwrite=True, nfiles=2)

        self.assertTrue(Path('tmp/Argo/{}/{}'.format('aoml', 3900407)).exists())

if __name__ == '__main__':
    Path('./tmp/Argo/aoml').mkdir(parents=True, exist_ok=True)
    Path('./tmp/NCEP/pres').mkdir(parents=True, exist_ok=True)
    Path('./tmp/NCEP/land').mkdir(parents=True, exist_ok=True)
    Path('./tmp/WOA18/o2sat').mkdir(parents=True, exist_ok=True)

    unittest.main()
    