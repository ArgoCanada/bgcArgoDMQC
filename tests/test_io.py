#!/usr/bin/python

import numpy as np
import pandas as pd
from netCDF4 import Dataset

from pathlib import Path
import unittest

import bgcArgoDMQC as bgc
class downloadTest(unittest.TestCase):

    def setUp(self):

        bgc.set_dirs(
            argo_path=Path(__file__).absolute().parent / 'test_data/Argo/dac',
            ncep_path=Path(__file__).absolute().parent / 'test_data/NCEP',
            woa_path=Path(__file__).absolute().parent / 'test_data/WOA18',
        )

        bgc.io.Path.ARGO_PATH.mkdir(exist_ok=True, parents=True)
        bgc.io.Path.NCEP_PATH.mkdir(exist_ok=True, parents=True)
        bgc.io.Path.WOA_PATH.mkdir(exist_ok=True, parents=True)


    def test_download_argo(self):

        wmo = 4901784
        bgc.io.get_argo(wmo, local_path=bgc.io.Path.ARGO_PATH, overwrite=False, nfiles=2)

        for mission in ['C', 'B', 'CB']:
            for mode in ['R', 'D', 'RD']:
                bgc.io.get_argo(wmo, local_path=bgc.io.Path.ARGO_PATH, overwrite=False, mission=mission, mode=mode, nfiles=1)
    
    def test_single_file_download(self):
        wmo = 4901784
        bgc.io.get_argo(f'/ifremer/argo/dac/meds/{wmo}/profiles/BD{wmo}_003.nc', local_path=bgc.io.Path.ARGO_PATH)

    def test_download_manipulate_argo(self):

        wmo = 4901784
        argo_path = Path(__file__).absolute().parent / 'test_data/Argo/dac'
        argo_path.mkdir(exist_ok=True)
        bgc.io.get_argo(wmo, local_path=argo_path, overwrite=True, nfiles=2)

        infile = list((argo_path / f'meds/{wmo}/profiles').glob('*.nc'))[0].as_posix()
        outfile = infile.replace('meds/', 'meds/E/')
        Path(outfile).parent.mkdir(parents=True, exist_ok=True)
        nc_out = bgc.io.iterate_dimension(infile, outfile, 'N_CALIB')
        nc = Dataset(infile, 'r')
        self.assertGreater(nc_out.dimensions['N_CALIB'].size, nc.dimensions['N_CALIB'].size)
        nc_out.close()

        new_var = dict(
            name='MY_NEW_VARIABLE',     # variable name, can be new or existing
            datatype=np.float64,        # variable datatype, can be np datatype or string
            dimensions=nc['PARAMETER_DATA_MODE'].dimensions,
            data=10*np.random.rand(*[nc.dimensions[d].size for d in nc['PARAMETER_DATA_MODE'].dimensions]),
            long_name='The new variable',
            standard_name='my_new_var',
            units='degree_celsius',
            valid_min=0,
            valid_max=1e9,
            resolution=0.001,
            _FillValue=99999.,
            comment='Added by Automated Test - not for Argo use'
        )

        nc.close()
        nc = bgc.io.append_variable(infile, new_var)

        self.assertTrue('MY_NEW_VARIABLE' in nc.variables.keys())
        nc.close()

    def test_profile_qc(self):

        test_flag_arrays = [
            pd.Series([4, 4, 4, 4]),
            pd.Series([1, 4, 4, 4, 4]),
            pd.Series([1, 4, 4, 4]),
            pd.Series([1, 1, 4, 4]),
            pd.Series([1, 1, 1, 4]),
            pd.Series([1, 1, 1, 1]),
        ]

        for a, g in zip(test_flag_arrays, ['F', 'E', 'D', 'C', 'B', 'A']):
            grade = bgc.io.profile_qc(a)
            self.assertEqual(grade, g)
    
    def test_index(self):
        
        core = bgc.io.read_index('core')
        traj = bgc.io.read_index('traj', source='local')
        traj = bgc.io.read_index('traj', source='remote')