#!/usr/bin/python

from pathlib import Path

import bgcArgoDMQC as bgc

import unittest

class downloadTest(unittest.TestCase):

    def test_download_refdata(self):

        bgc.io.get_ncep('rhum', local_path='tmp/NCEP', overwrite=True, years=[2019, 2020])
        bgc.io.get_woa18('O2sat', local_path='tmp/WOA18', nfiles=0)

    def test_download_argo(self):

        wmo = 4901784
        bgc.io.get_argo(wmo, local_path='tmp/Argo/dac/', overwrite=True, nfiles=2)

if __name__ == '__main__':

    unittest.main()
    