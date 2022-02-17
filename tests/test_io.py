import smps
import pandas as pd
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import pytest
from smps.io import load_sample
import unittest

class TestClass(unittest.TestCase):

    def test_generic(self):
        res = smps.io.smps_from_txt("https://raw.githubusercontent.com/dhhagan/py-smps/master/sample-data/boston_wintertime.txt", column=False)

        # create a model
        m = smps.models.GenericParticleSizer(data=res['data'], bins=res['bins'],
                                                dp_units='nm', fmt='dndlogdp')

        assert isinstance(m, smps.models.GenericParticleSizer)
        assert m.bins.shape[0] == len(m.s_multiplier)
        assert m.bins.shape[0] == len(m.v_multiplier)
        assert m.bins.shape[0] == len(m.dlogdp)

        # test .copy()
        cpy = m.copy()

        assert isinstance(cpy, smps.models.GenericParticleSizer)
        assert cpy != m

        # test the scan stats
        stats = m.scan_stats
        assert isinstance(stats, pd.DataFrame)
        assert list(stats.columns) == (list(set(m.data.columns) - set(m.dn.columns)))

        # test stats
        stats = m.stats(weight='number')

        self.assertIsInstance(stats, pd.DataFrame)
        self.assertTrue('number' in stats.columns)
        self.assertTrue('surface_area' in stats.columns)
        self.assertTrue('volume' in stats.columns)
        self.assertTrue('mass' in stats.columns)
        self.assertTrue('AM' in stats.columns)
        self.assertTrue('GM' in stats.columns)
        self.assertTrue('GSD' in stats.columns)
        self.assertTrue('Mode' in stats.columns)

        pm1 = m.integrate(weight='number', dmin=0, dmax=1.)
        cut1 = m.integrate(weight='number', dmin=0, dmax=0.5)

        self.assertGreaterEqual(pm1.sum(), cut1.sum())

        pm1 = m.integrate(weight='surface', dmin=0, dmax=1.)
        cut1 = m.integrate(weight='surface', dmin=0, dmax=0.5)

        self.assertGreaterEqual(pm1.sum(), cut1.sum())

        pm1 = m.integrate(weight='volume', dmin=0, dmax=1.)
        cut1 = m.integrate(weight='volume', dmin=0, dmax=0.5)

        self.assertGreaterEqual(pm1.sum(), cut1.sum())

        pm1 = m.integrate(weight='mass', dmin=0, dmax=1., rho=1.)
        pm2 = m.integrate(weight='mass', dmin=0, dmax=1., rho=1.65)

        self.assertGreaterEqual(pm2.sum(), pm1.sum())

        sliced = m.slice(start=(m.data.index[0] + pd.Timedelta(hours=12)), inplace=False)

        self.assertLessEqual(sliced.data.shape[0], m.data.shape[0])

        rs = m.resample("30min", inplace=False)

        self.assertGreaterEqual(m.data.shape[0], rs.data.shape[0])

    def test_generic_calculations(self):
        number = smps.io.smps_from_txt(
            "https://raw.githubusercontent.com/dhhagan/py-smps/master/tests/datafiles/test_data_number.txt",
            column=False, 
            as_dict=False
        )

        surface = smps.io.smps_from_txt(
            "https://raw.githubusercontent.com/dhhagan/py-smps/master/tests/datafiles/test_data_surface_area.txt",
            column=False, 
            as_dict=False
        )

        volume = smps.io.smps_from_txt(
            "https://raw.githubusercontent.com/dhhagan/py-smps/master/tests/datafiles/test_data_volume.txt",
            column=False, 
            as_dict=False
        )

        # using the number data, compare our calculations using 'stats()' to
        # those that the AIM software output
        stats = number.stats()

        cols_to_check = []
        cols_to_check.append(("number", "Total Conc.(#/cm³)"))
        cols_to_check.append(("AM", "Mean(nm)"))
        cols_to_check.append(("GM", "Geo. Mean(nm)"))
        cols_to_check.append(("Mode", "Mode(nm)"))
        cols_to_check.append(("GSD", "Geo. Std. Dev."))

        print ("")

        for pkg_col, aim_col in cols_to_check:
            m, b, r, p, e = linregress(
                                stats[pkg_col].values,
                                number.scan_stats[aim_col].values)
            
            print ("\t{} v {}: \n\t\tr2={:.3f}, m={:.2f}".format(pkg_col, aim_col, r**2, m))

            # make sure the correlation is above 0.99
            self.assertGreaterEqual(r**2, 0.99)

            # make sure the slope is between 0.99 and 1.01
            self.assertGreaterEqual(m, 0.99)
            self.assertLessEqual(m, 1.01)

        # check the surface-area-weighted stats
        stats = number.stats(weight='surface')

        cols_to_check = []
        cols_to_check.append(("surface_area", "Total Conc.(nm²/cm³)"))
        cols_to_check.append(("AM", "Mean(nm)"))
        cols_to_check.append(("GM", "Geo. Mean(nm)"))
        cols_to_check.append(("Mode", "Mode(nm)"))
        cols_to_check.append(("GSD", "Geo. Std. Dev."))

        for pkg_col, aim_col in cols_to_check:
            if pkg_col in ['surface_area']:
                stats[pkg_col]*=1e6

            m, b, r, p, e = linregress(
                                stats[pkg_col].values,
                                surface.scan_stats[aim_col].values)

            print("\t{} v {}: \n\t\tr2={:.3f}, m={:.2f}".format(pkg_col, aim_col, r**2, m))

            # make sure the correlation is above 0.99
            self.assertGreaterEqual(r**2, 0.975)

            # make sure the slope is between 0.99 and 1.01
            self.assertGreaterEqual(m, 0.99)
            self.assertLessEqual(m, 1.01)


        # check the volume-weighted stats
        stats = number.stats(weight='volume')

        cols_to_check = []
        cols_to_check.append(("volume", "Total Conc.(nm³/cm³)"))
        cols_to_check.append(("AM", "Mean(nm)"))
        cols_to_check.append(("GM", "Geo. Mean(nm)"))
        cols_to_check.append(("Mode", "Mode(nm)"))
        cols_to_check.append(("GSD", "Geo. Std. Dev."))

        for pkg_col, aim_col in cols_to_check:
            if pkg_col in ['volume']:
                stats[pkg_col]*=1e9

            m, b, r, p, e = linregress(
                                stats[pkg_col].values,
                                volume.scan_stats[aim_col].values)

            print("\t{} v {}: \n\t\tr2={:.3f}, m={:.2f}".format(pkg_col, aim_col, r**2, m))
            # make sure the correlation is above 0.99
            self.assertGreaterEqual(r**2, 0.975)

            # make sure the slope is between 0.99 and 1.01
            self.assertGreaterEqual(m, 0.95)
            self.assertLessEqual(m, 1.05)

    def test_models_smps(self):
        pass

#    def test_models_alphasense_opcn2(self):
#        opcn2 = smps.io.opcn2_from_text(
#                    os.path.join(datadir, "OPC2_001.CSV"),
#                    as_dict=False)
#
#        self.assertIsInstance(opcn2, smps.models.AlphasenseOpcN2)
#
#        # test the mass calculation - should match what Alphasense gets
#        cols_to_check = []
#        cols_to_check.append(("PM1(ug/m3)", 1))
#        cols_to_check.append(("PM2.5", 2.5))
#        cols_to_check.append(("PM10", 10))
#
#        print ()
#        for opc_col, dmax in cols_to_check:
#            m, b, r, p, e = linregress(
#                                opcn2.scan_stats[opc_col].values,
#                                opcn2.integrate(weight='mass', rho=1.65, dmax=dmax).values)
#
#            print ("\t{} @ dmax={}".format(opc_col, dmax))
#            print ("\t\tm={:.3f}, r2={:.3f}".format(m, r**2))

    def test_models_modulairpm(self):
        path_to_modpm = "https://raw.githubusercontent.com/quant-aq/py-smps/master/tests/datafiles/MOD-PM-SAMPLE.csv"
        modpm_df = pd.read_csv(path_to_modpm)

        modpm_model = smps.models.ModulairPM(data=modpm_df, fmt='dlogdn')
        
        integrated = modpm_model.integrate(weight='number', dmin=0.05, dmax=0.5)

    def test_models_modulair(self):
        path_to_mod = "https://raw.githubusercontent.com/quant-aq/py-smps/master/tests/datafiles/MODULAIR-SAMPLE.csv"
        mod_df = pd.read_csv(path_to_mod)

        mod_model = smps.models.Modulair(data=mod_df)

        integrated = mod_model.integrate(weight='number', dmin=0.05, dmax=0.5)
    
    def test_models_pops(self):
        pass

    def test_fit(self):
        from smps.fit import LogNormal

        r = smps.io.smps_from_txt(
            "https://raw.githubusercontent.com/dhhagan/py-smps/master/tests/datafiles/test_data_number.txt",
            column=False, 
            as_dict=False
        )

        # fit 1 mode in volume number space
        m = LogNormal()

        results = m.fit(X=r.midpoints, Y=r.dndlogdp.mean(), modes=1)

    def test_calculations_with_nans(self):
        number = smps.io.smps_from_txt(
            "https://raw.githubusercontent.com/dhhagan/py-smps/master/tests/datafiles/test_data_number.txt",
            column=False, 
            as_dict=False
        )

        number.resample("1min", inplace=True)

        stats = number.stats(weight='number')

def testObjects():
    obj = load_sample("boston")

    assert isinstance(obj, smps.SMPS)
    assert isinstance(obj, smps.GenericParticleSizer)
    assert isinstance(obj.data, pd.DataFrame)

    # Check the keys
    assert hasattr(obj, "meta")
    assert isinstance(obj.meta, dict)
    assert isinstance(obj.meta.get('units'), str)
    assert isinstance(obj.meta.get('weight'), str)
    assert isinstance(obj.data, pd.DataFrame)
    assert isinstance(obj.bins, np.ndarray)
    assert isinstance(obj.bin_labels, list)
