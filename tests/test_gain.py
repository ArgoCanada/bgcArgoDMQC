#!/usr/bin/python

import bgcArgoDMQC as bgc
from pathlib import Path
import unittest

class gainTest(unittest.TestCase):

    def test_gain(self):

        ncep_path = Path('test_data/NCEP').absolute()
        woa_path = Path('test_data/WOA18').absolute()

        ncep_path.mkdir(exist_ok=True)
        woa_path.mkdir(exist_ok=True)

        bgc.io.get_ncep('rhum', local_path=ncep_path, overwrite=True, years=[2019, 2020])
        bgc.io.get_ncep('pres', local_path=ncep_path, overwrite=True, years=[2019, 2020])
        bgc.io.get_ncep('land', local_path=ncep_path, overwrite=True, years=[2019, 2020])
        bgc.io.get_woa18('O2sat', local_path=woa_path)

        wmo = 6902807
        argo_path = Path('test_data/Argo/dac').absolute()
        argo_path.mkdir(exist_ok=True)
        bgc.io.get_argo(wmo, local_path=argo_path, overwrite=True, nfiles=5)

        bgc.set_dirs(
            argo_path=argo_path,
            woa_path=woa_path,
            ncep_path=ncep_path,
        )

        flt = bgc.sprof(wmo)
        woa_gains = flt.calc_gains(ref='WOA')
        ncep_gains = flt.calc_gains()

        flt.plot('gains', ref='WOA')
        flt.plot('gains', ref='NCEP')

if __name__ == '__main__':

    unittest.main()