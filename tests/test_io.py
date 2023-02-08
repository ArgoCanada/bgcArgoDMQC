#!/usr/bin/python

import bgcArgoDMQC as bgc
from pathlib import Path
import unittest

class downloadTest(unittest.TestCase):

    def test_download_refdata(self):

        ncep_path = Path(__file__).absolute().parent / 'test_data/NCEP'
        woa_path = Path(__file__).absolute().parent / 'test_data/WOA18'

        ncep_path.mkdir(exist_ok=True)
        woa_path.mkdir(exist_ok=True)

        bgc.io.get_ncep('rhum', local_path=ncep_path, overwrite=True, years=[2019, 2020])
        bgc.io.get_woa18('O2sat', local_path=woa_path, nfiles=0)

    def test_download_manipulate_argo(self):

        wmo = 4901784
        argo_path = Path(__file__).absolute().parent / 'test_data/Argo/dac'
        argo_path.mkdir(exist_ok=True)
        bgc.io.get_argo(wmo, local_path=argo_path, overwrite=True, nfiles=2)

        infile = list((argo_path / f'meds/{wmo}/profiles').glob('*.nc'))[-1].as_posix()
        outfile = infile.replace('.nc', '_out.nc')
        bgc.io.iterate_dimension(infile, outfile, 'N_CALIB')
