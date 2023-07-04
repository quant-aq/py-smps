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
        
    # def test_variable_kappa(self):
    #     path = "https://raw.githubusercontent.com/quant-aq/py-smps/master/tests/datafiles/MOD-PM-SAMPLE.csv"
    #     df = pd.read_csv(path)

    #     # Create the object
    #     obj = smps.models.ModulairPM(data=df)
        
    #     # First, compute PM with a constant density
    #     pm = obj.integrate(weight='mass', dmin=0., dmax=1., kappa=0.3)
        
    #     # Compute PM with a size-dependant density
    #     def custom_kappa(dp):
    #         if dp < 2.5:
    #             return .3
            
    #         return 0.
        
    #     pm = obj.integrate(weight='mass', dmin=0., dmax=10., kappa=custom_kappa)
