import unittest
import smps
import os
import pandas as pd
import numpy as np

basedir = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(basedir, "datafiles")

class SetupTestCase(unittest.TestCase):
    def setUp(self):
        pass

        # Import some test data
        #self.data_number = load_file(os.path.join(basedir, "datafiles/test_data_number.txt"), column=False)
        #self.data_diameter = load_file(os.path.join(basedir, "datafiles/test_data_diameter.txt"), column=False)
        #self.data_surface_area = load_file(os.path.join(basedir, "datafiles/test_data_surface_area.txt"), column=False)
        #self.data_volume = load_file(os.path.join(basedir, "datafiles/test_data_volume.txt"), column=False)

    def tearDown(self):
        pass

    def test_smps_from_txt(self):
        res = smps.io.smps_from_txt(
                        os.path.join(datadir, "test_data_number.txt"),
                        column=False)

        self.assertTrue('meta' in res.keys())
        self.assertTrue('units' in res.keys())
        self.assertTrue('weight' in res.keys())
        self.assertTrue('data' in res.keys())
        self.assertTrue('bins' in res.keys())
        self.assertTrue('bin_labels' in res.keys())
        self.assertTrue('bin_prefix' in res.keys())

        self.assertTrue(isinstance(res['meta'], dict))
        self.assertTrue(isinstance(res['units'], str))
        self.assertTrue(isinstance(res['weight'], str))
        self.assertTrue(isinstance(res['data'], pd.DataFrame))
        self.assertTrue(isinstance(res['bins'], np.ndarray))
        self.assertTrue(isinstance(res['bin_labels'], list))
        self.assertTrue(isinstance(res['bin_prefix'], str))

    def test_smps_from_txt_columns(self):
        res = smps.io.smps_from_txt(
                        os.path.join(datadir, "test_data_column.txt"),
                        column=True)

        self.assertTrue('meta' in res.keys())
        self.assertTrue('units' in res.keys())
        self.assertTrue('weight' in res.keys())
        self.assertTrue('data' in res.keys())
        self.assertTrue('bins' in res.keys())
        self.assertTrue('bin_labels' in res.keys())
        self.assertTrue('bin_prefix' in res.keys())

        self.assertTrue(isinstance(res['meta'], dict))
        self.assertTrue(isinstance(res['units'], str))
        self.assertTrue(isinstance(res['weight'], str))
        self.assertTrue(isinstance(res['data'], pd.DataFrame))
        self.assertTrue(isinstance(res['bins'], np.ndarray))
        self.assertTrue(isinstance(res['bin_labels'], list))
        self.assertTrue(isinstance(res['bin_prefix'], str))

    def test_models_generic_methods(self):
        res = smps.io.smps_from_txt(
                        os.path.join(datadir, "test_data_number.txt"),
                        column=False)

        # create a model
        m = smps.models.GenericParticleSizer(
                            data=res['data'],
                            bins=res['bins'],
                            dp_units='nm',
                            fmt='dndlogdp')

        # make sure the model has all of the methods!
        self.assertIsInstance(m, smps.models.GenericParticleSizer)

        self.assertEqual(m.bins.shape[0], len(m.s_multiplier))
        self.assertEqual(m.bins.shape[0], len(m.v_multiplier))
        self.assertEqual(m.bins.shape[0], len(m.dlogdp))

        # test .copy()
        cpy = m.copy()

        self.assertIsInstance(cpy, smps.models.GenericParticleSizer)
        self.assertNotEqual(cpy, m)

        # test scan stats
        stats = m.scan_stats
        self.assertIsInstance(stats, pd.DataFrame)
        self.assertEqual(list(stats.columns), list(set(m.data.columns) - set(m.dn.columns)))

        # test ._subselect_frame()
        

    def test_models_smps(self):
        pass

    def test_models_alphasense_opcn2(self):
        pass

    def test_generic_calculations(self):
        pass

"""
    def test_smps_copy(self):
        df = self.data_number

        cpy = df.copy()

        self.assertTrue(df.raw.equals(cpy.raw))

    def test_smps_calculations(self):
        s = self.data_number

        # Test the first bin of dlogdp
        self.assertEqual(len(s.dlogdp), s.bins.shape[0])
        self.assertTrue(s.dndlogdp.equals(s.raw[s.bin_labels]))

        # Make sure that dDdlogDp is correct
        _calculated = self.data_number.dddlogdp * 1e-3 # divide by 1000 to go from um to nm
        _reference = self.data_diameter.dndlogdp # nanometers

        self.assertEqual(round(_calculated.iloc[0][0], 2), round(_reference.iloc[0][0], 2))

        # Make sure that dSdlogDp is correct

    def test_datatypes(self):
        df = self.data_number

        self.assertEqual(df.raw['Median'].dtype, float)
        self.assertEqual(df.raw['Mean'].dtype, float)
        self.assertEqual(df.raw['Mode'].dtype, float)
        self.assertEqual(df.raw['GM'].dtype, float)
        self.assertEqual(df.raw['GSD'].dtype, float)
        self.assertEqual(df.raw['Total Conc.'].dtype, float)

    def test_resampling(self):
        df = self.data_number

        #self.assertIsNotNone(df.resample('5min', inplace=False))
        self.assertTrue(df.resample('5min', inplace=True))

    def test_smps_model(self):
        pass
        #model = load_sample('boston')

        # Check dlo

    def test_stats(self):
        df = self.data_number
        df2 = self.data_surface_area
        df3 = self.data_volume

        # Retrieve the stats
        stats = df.stats(weight='number')

        self.assertTrue('Total Number' in stats.columns)
        self.assertTrue('Total Surface Area' in stats.columns)
        self.assertTrue('Total Volume' in stats.columns)
        self.assertTrue('Mean' in stats.columns)

        # Make sure the GM, GSD, and Mean are all within 1% error
        def one_pct_error(x1, x2):
            diff = abs(x1 - x2) / x1

            return True if diff <= 0.001 else False

        self.assertTrue(one_pct_error(stats["GM"][0], df.scan_stats['GM'][0]))
        self.assertTrue(one_pct_error(stats["Mean"][0], df.scan_stats['Mean'][0]))
        self.assertTrue(one_pct_error(stats["GSD"][0], df.scan_stats['GSD'][0]))

        # Repeat for Surface-Area weighted Statistics
        stats = df.stats(weight='surface_area')

        self.assertTrue(one_pct_error(stats["GM"][0], df2.scan_stats['GM'][0]))
        self.assertTrue(one_pct_error(stats["Mean"][0], df2.scan_stats['Mean'][0]))
        self.assertTrue(one_pct_error(stats["GSD"][0], df2.scan_stats['GSD'][0]))

        # Repeat for Volume weighted Statistics
        stats = df.stats(weight='volume')

        self.assertTrue(one_pct_error(stats["GM"][0], df3.scan_stats['GM'][0]))
        self.assertTrue(one_pct_error(stats["Mean"][0], df3.scan_stats['Mean'][0]))
        self.assertTrue(one_pct_error(stats["GSD"][0], df3.scan_stats['GSD'][0]))

    def test_import_aim10_2(self):
        df = load_file(os.path.join(basedir, "datafiles/aim10_2.txt"), column=False)
"""
