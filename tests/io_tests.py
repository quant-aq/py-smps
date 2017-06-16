import unittest
import smps
from smps.io import load_sample

class SetupTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_load_sample(self):
        # Load the sample
        df = load_sample('boston')

        self.assertIsInstance(df, smps.io.SMPS)

    def test_load_sample_column(self):
        df = load_sample('chamber')

        self.assertIsInstance(df, smps.io.SMPS)

    def test_datatypes(self):
        df = load_sample('boston')

        self.assertEqual(df.raw['Median'].dtype, float)
        self.assertEqual(df.raw['Mean'].dtype, float)
        self.assertEqual(df.raw['Mode'].dtype, float)
        self.assertEqual(df.raw['GM'].dtype, float)
        self.assertEqual(df.raw['GSD'].dtype, float)
        self.assertEqual(df.raw['Total Conc.'].dtype, float)

    def test_resampling(self):
        df = load_sample('boston')

        rs = df.raw.resample('5min').mean()

        self.assertIsNotNone(rs['Mean'])


    def test_smps_model(self):
        model = load_sample('boston')

        # Check dlo
