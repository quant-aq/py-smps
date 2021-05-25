import smps
import os
import pandas as pd
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import pytest
from smps.io import load_sample


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
