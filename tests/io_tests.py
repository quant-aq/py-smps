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

    def test_datatypes(self):
        df = load_sample('boston')

        self.assertEqual(df.raw['Median'].dtype, float)

    def test_smps_model(self):
        model = load_sample('boston')

        # Check dlo
