#!/usr/bin/python

import bgcArgoDMQC as bgc
from pathlib import Path
import unittest

class downloadTest(unittest.TestCase):

    def test_download_refdata(self):

        ncep_path = Path('test_data/NCEP').absolute()
        woa_path = Path('test_data/WOA18').absolute()

        ncep_path.mkdir(exist_ok=True)
        woa_path.mkdir(exist_ok=True)

        bgc.io.get_ncep('rhum', local_path=ncep_path, overwrite=True, years=[2019, 2020])
        bgc.io.get_woa18('O2sat', local_path=woa_path, nfiles=0)

    def test_download_argo(self):

        wmo = 4901784
        argo_path = Path('test_data/Argo/dac').absolute()
        argo_path.mkdir(exist_ok=True)
        bgc.io.get_argo(wmo, local_path=argo_path, overwrite=True, nfiles=2)

if __name__ == '__main__':

    unittest.main()