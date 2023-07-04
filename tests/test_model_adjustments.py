import smps
import pandas as pd
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import pytest
from smps.io import load_sample
import unittest

class TestClass(unittest.TestCase):

    def test_variable_density(self):
        path = "https://raw.githubusercontent.com/quant-aq/py-smps/master/tests/datafiles/MOD-PM-SAMPLE.csv"
        df = pd.read_csv(path)

        # Create the object
        obj = smps.models.ModulairPM(data=df)
        
        # First, compute PM with a constant density
        pm = obj.integrate(weight='mass', dmin=0., dmax=1., rho=1.)
        
        # Compute PM with a size-dependant density
        def custom_density(dp):
            if dp < 2.5:
                return 1.
            
            return 2.
        
        pm = obj.integrate(weight='mass', dmin=0., dmax=10., rho=custom_density)
        
    def test_variable_kappa(self):
        path = "https://raw.githubusercontent.com/quant-aq/py-smps/master/tests/datafiles/MOD-PM-SAMPLE.csv"
        df = pd.read_csv(path, parse_dates=['timestamp']).set_index("timestamp")

        # Create the object
        obj = smps.models.ModulairPM(data=df)
        
        # If kappa is set, but rh is not, raise an exception
        with self.assertRaises(AttributeError):
            pm = obj.integrate(weight='mass', dmin=0., dmax=1., kappa=0.3)
            
        # If the rh column doesn't exist, raise an exception
        with self.assertRaises(smps.models.ValidationError):
            pm = obj.integrate(weight='mass', dmin=0., dmax=1., kappa=0.3, rh='rh')
            
        # Compute with kappa and a real rh (mass)
        pm = obj.integrate(weight='mass', dmin=0., dmax=2.5, kappa=0.3, rh="sample_rh")
        
        # Compute with kappa and a real rh (number)
        pm = obj.integrate(weight='number', dmin=0., dmax=2.5, kappa=0.3, rh="sample_rh")
        
        # Compute with kappa and a real rh (surface)
        pm = obj.integrate(weight='surface', dmin=0., dmax=2.5, kappa=0.3, rh="sample_rh")
        
        # Compute with kappa and a real rh (volume)
        pm = obj.integrate(weight='volume', dmin=0., dmax=2.5, kappa=0.3, rh="sample_rh")
        
        # make sure the values are lower at a higher kappa
        pmx = obj.integrate(weight='mass', dmin=0., dmax=2.5, kappa=0.3, rh="sample_rh")
        pmy = obj.integrate(weight='mass', dmin=0., dmax=2.5, kappa=1, rh="sample_rh")
        
        self.assertLess(pmy.mean(), pmx.mean())