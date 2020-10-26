#!/usr/bin/python

from pathlib import Path
import shutil

from netCDF4 import Dataset

import matplotlib.pyplot as plt

import unittest
import numpy as np
import pandas as pd
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

        bgc.io.get_ncep('pres', local_path='tmp/NCEP', overwrite=True, years=[2019, 2020])
        bgc.io.get_ncep('land', local_path='tmp/NCEP', overwrite=True)
        bgc.io.get_ncep('rhum', local_path='tmp/NCEP', overwrite=True, years=[2019, 2020])

        self.assertTrue(Path('tmp/NCEP/pres/pres.sfc.gauss.2019.nc').exists())
        self.assertTrue(Path('tmp/NCEP/land/land.sfc.gauss.nc').exists())
        self.assertTrue(Path('tmp/NCEP/rhum/rhum.sig995.2019.nc').exists())

    def test_download_woa(self):

        bgc.io.get_woa18('O2sat', local_path='tmp/WOA18', overwrite=True)

        self.assertTrue(Path('tmp/WOA18/o2sat/woa18_all_O00_01.nc').exists())

    def test_download_argo(self):

        dac = [bgc.io.get_dac(w) for w in [wmo-1,wmo]]

        dacpath = '/ifremer/argo/dac'
        fltpath = ['{}/{}/{}'.format(dacpath, d, w) for d, w in zip(dac, [wmo-1,wmo])]
        bgc.io.get_argo(fltpath, local_path='tmp/Argo', overwrite=True)
        bgc.io.get_argo('/ifremer/argo/dac/aoml', [3900407, 4900345], local_path='tmp/Argo/aoml', overwrite=True)

        self.assertTrue(Path('tmp/Argo/{}/{}'.format('meds',wmo-1)).exists())
        self.assertTrue(Path('tmp/Argo/{}/{}'.format('meds',wmo)).exists())
        self.assertTrue(Path('tmp/Argo/{}/{}'.format('aoml',3900407)).exists())
        self.assertTrue(Path('tmp/Argo/{}/{}'.format('aoml',4900345)).exists())

class sprofTest(unittest.TestCase):

    def test_sprof(self):
        sprof = bgc.sprof(wmo)
        sprof.clean()
        sprof.reset()
        sprof.clean(bad_flags=[3,4])
        df = sprof.to_dataframe()
        sprof.reset()

        self.assertIsInstance(sprof, bgc.sprof)
        self.assertIs(type(df), pd.core.frame.DataFrame)

    def test_calc_gains(self):
        sprof = bgc.sprof(wmo)
        ncep_gains = sprof.calc_gains()
        woa_gains  = sprof.calc_gains(ref='WOA')

        self.assertIs(type(ncep_gains), np.ndarray)
        self.assertIs(type(woa_gains), np.ndarray)

class profilesTest(unittest.TestCase):

    def test_profiles(self):
        prof = bgc.profiles(wmo)
        prof.clean()
        prof.reset()
        prof.clean(bad_flags=4)
        prof.reset()
        # test multiple profs
        profs = bgc.profiles([wmo-1, wmo])
        # test specific cycles
        cycs = bgc.profiles(wmo, cycles=np.arange(1,10))

        df = cycs.to_dataframe()

        self.assertIsInstance(prof, bgc.profiles)
        self.assertIsInstance(profs, bgc.profiles)
        self.assertIsInstance(cycs, bgc.profiles)
        self.assertIs(type(df), pd.core.frame.DataFrame)

    def test_calc_gains(self):
        prof = bgc.profiles(wmo)
        woa_gains  = prof.calc_gains(ref='WOA')

        self.assertIs(type(woa_gains), np.ndarray)

class plottingTest(unittest.TestCase):

    def test_gain_plot(self):

        syn  = bgc.sprof(4902480)
        prof = bgc.profiles([4902480, 6902905])

        syn.calc_gains()
        syn.calc_gains(ref='WOA')

        g_ncep = syn.plot('gain', ref='NCEP')
        self.assertIsInstance(g_ncep, bgc.fplt.pltClass)
        plt.close(g_ncep.fig)

        g_woa  = syn.plot('gain', ref='WOA')
        self.assertIsInstance(g_woa, bgc.fplt.pltClass)
        plt.close(g_woa.fig)

    def test_scatter_plot(self):

        syn = bgc.sprof(4902480)
        g = syn.plot('cscatter', varname='DOXY', ylim=(0,500))

        self.assertIsInstance(g, bgc.fplt.pltClass)

        plt.close(g.fig)

    def test_profile_plot(self):
        syn = bgc.sprof(4902480)

        g_pres = syn.plot('profiles', varlist=['TEMP', 'DOXY'])
        self.assertIsInstance(g_pres, bgc.fplt.pltClass)
        plt.close(g_pres.fig)

        g_pden = syn.plot('profiles', varlist=['PSAL', 'DOXY'], Nprof=5, Ncycle=3, zvar='PDEN')
        self.assertIsInstance(g_pden, bgc.fplt.pltClass)
        plt.close(g_pden.fig)

class coreTest(unittest.TestCase):

    bgc.io.check_index(mode='install')
    bgc.io.check_index()

    bgc.io.update_index()

    def test_index_files(self):
        # get index files
        bgc_index  = bgc.get_index('bgc')
        core_index = bgc.get_index('global')
        syn_index  = bgc.get_index('synthetic')

        self.assertIs(type(bgc_index), pd.core.frame.DataFrame)
        self.assertIs(type(core_index), pd.core.frame.DataFrame)
        self.assertIs(type(syn_index), pd.core.frame.DataFrame)

    def test_qc_read(self):
        # read QC test
        nc = Dataset(Path('tmp/Argo/aoml/3900407/profiles/D3900407_188.nc'))
        qcp, qcf = bgc.read_history_qctest(nc)
        bgc.util.display_qctests(qcp, qcf)

        self.assertIs(type(qcp), np.str_)
        self.assertIs(type(qcf), np.str_)

    def test_information_criteria(self):
        # aic and bic
        data  = np.random.randn(30)
        resid = np.random.randn(20)
        aic = bgc.util.bic(data, resid)

        self.assertIs(type(aic), np.float)

        data  = np.random.randn(30)
        resid = np.random.randn(20)
        bic = bgc.util.bic(data, resid)

        self.assertIs(type(bic), np.float)   

    def test_response_time_correction(self):
        # time correction
        time = np.arange(0,20,1)
        temp = 12*np.random.rand(time.shape[0])
        doxy = 200*np.random.rand(time.shape[0])
        thickness = 200

        doxy_adj = bgc.correct_response_time(time, doxy, temp, 120)
        self.assertIs(type(doxy_adj), np.ndarray)

        doxy_adj_Tconst = bgc.correct_response_time_Tconst(time, doxy, 70)
        self.assertIs(type(doxy_adj_Tconst), np.ndarray)

    def test_get_var_by(self):
        # get var by
        dv = dict(a=np.random.randn(20), b=np.random.randn(20))
        v = bgc.util.get_var_by('a', 'b', dv)

        self.assertIs(type(v), np.ndarray)

class otherTest(unittest.TestCase):

    def test_pO2(self):
        syn = bgc.sprof(4902480)
        pO2 = bgc.unit.pO2(syn.DOXY, syn.PSAL, syn.TEMP)
        self.assertIs(type(pO2), np.ndarray)

    def test_read_gain_value(self):
        g, eq, comment = bgc.util.read_gain_value(Path('tmp/Argo/aoml/4900345/profiles/BD4900345_024.nc'))

        self.assertIs(type(g[0]), np.str_)
        self.assertIs(type(eq[0]), np.str_)
        self.assertIs(type(comment[0]), np.str_)

    def test_unit_conversion(self):
        S = np.random.rand(20)
        T = np.random.rand(20)

        sol = bgc.unit.oxy_sol(S, T, unit='millimole/m3')

        self.assertIs(type(sol), np.ndarray)

        with self.assertRaises(ValueError):
            sol = bgc.unit.oxy_sol(S, T, unit='deg C')

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
    
    shutil.rmtree('tmp')