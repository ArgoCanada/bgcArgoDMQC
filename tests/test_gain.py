#!/usr/bin/python

import bgcArgoDMQC as bgc
from pathlib import Path
import unittest

class gainTest(unittest.TestCase):

    def test_gain(self):

        ncep_path = Path(__file__).absolute().parent / 'test_data/NCEP'
        woa_path = Path(__file__).absolute().parent / 'test_data/WOA18'

        ncep_path.mkdir(exist_ok=True)
        woa_path.mkdir(exist_ok=True)

        bgc.io.get_ncep('rhum', local_path=ncep_path, overwrite=True, years=[2019, 2022])
        bgc.io.get_ncep('pres', local_path=ncep_path, overwrite=True, years=[2019, 2022])
        bgc.io.get_ncep('land', local_path=ncep_path, overwrite=True, years=[2019, 2022])
        bgc.io.get_woa18('O2sat', local_path=woa_path)

        wmo = 6902870
        argo_path = Path(__file__).absolute().parent / 'test_data/Argo/dac'
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

        flt.plot('gain', ref='WOA')
        flt.plot('gain', ref='NCEP')

        flt.update_field('DOXY_ADJUSTED', flt.gain*flt.DOXY)
        flt.update_field('DOXY_ADJUSTED_QC', flt.__fillvalue__['DOXY'], where=flt.DOXY_ADJUSTED_QC == 4)
        flt.export_files()
