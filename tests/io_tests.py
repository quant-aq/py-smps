import unittest
import smps
from smps.io import load_sample, load_file
import os

basedir = os.path.dirname(os.path.abspath(__file__))

class SetupTestCase(unittest.TestCase):
    def setUp(self):

        # Import some test data
        self.data_number = load_file(os.path.join(basedir, "datafiles/test_data_number.txt"), column=False)
        self.data_diameter = load_file(os.path.join(basedir, "datafiles/test_data_diameter.txt"), column=False)
        self.data_surface_area = load_file(os.path.join(basedir, "datafiles/test_data_surface_area.txt"), column=False)
        self.data_volume = load_file(os.path.join(basedir, "datafiles/test_data_volume.txt"), column=False)

    def tearDown(self):
        pass

    def test_load_sample(self):
        # Load the sample
        df = load_sample('boston')

        self.assertIsInstance(df, smps.io.SMPS)

    def test_load_sample_column(self):
        df = load_sample('chamber')

        self.assertIsInstance(df, smps.io.SMPS)

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
