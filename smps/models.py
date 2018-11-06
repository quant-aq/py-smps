#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import math
import copy
import joblib
import json

class GenericParticleSizer(object):
    """
        GenericParticleSizer is the base class for all
        available particle sizing instruments. It contains
        all of the basic functionality and methods used for
        making calculations.

        It can receive data in several formats including: dN, dNdlogDp
    """
    def __init__(self, data, bins, **kwargs):
        """
        :param data:
        :param bins:
        :param dp_units: ['um', 'nm']
        :param fmt: ['dn', 'dndlogdp']
        """
        self.data = data
        self.bins = bins
        self.meta = kwargs.pop('meta', dict())
        self.bin_labels = kwargs.pop('bin_labels', None)

        self.meta['weight'] = kwargs.pop("weight", "number")
        self.meta['units'] = kwargs.pop("units", "dw/dlogdp")

        # if bin_labels were not provided, try to generate them
        if not self.bin_labels:
            self.bin_labels = [c for c in self.data.columns if kwargs.pop("bin_prefix", "bin") in c]

        # if no bin labels exist, raise an error
        assert(len(self.bin_labels) > 0), "No bin labels have been found or specified."

        # make sure bins is has three columns
        assert(self.bins.shape[1] == 3), "`bins` must be an Nx3 array."

        # we always work with our units in um, so convert
        if kwargs.pop('dp_units', 'um') == 'nm':
            self.bins = self.bins * 1e-3

        # if the data is in a format other than [weight='number', units='dw/dlogdp']
        if kwargs.pop('fmt', 'dn') == 'dn':
            self.data[self.bin_labels] = self.data[self.bin_labels].div(self.dlogdp)

    @property
    def s_multiplier(self):
        return 4 * np.pi * (self.bins[:, 1]/2)**2

    @property
    def v_multiplier(self):
        return (4./3)*np.pi*(self.bins[:, 1]/2)**3

    @property
    def dlogdp(self):
        """"""
        return np.log10(self.bins[:, -1]) - np.log10(self.bins[:, 0])

    @property
    def dn(self):
        """"""
        return self.dndlogdp.mul(self.dlogdp)

    @property
    def dndlogdp(self):
        """Return the histogram in units of dNd/logDp [#/cm3]"""
        return self.data[self.bin_labels]

    @property
    def dd(self):
        """"""
        return self.dddlogdp.mul(self.dlogdp)

    @property
    def dddlogdp(self):
        """Return the histogram in units of dDpd/logDp [um/cm3]"""
        return self.data[self.bin_labels].mul(self.midpoints)

    @property
    def ds(self):
        """"""
        return self.dsdlogdp.mul(self.dlogdp)

    @property
    def dsdlogdp(self):
        """Return the histogram in units of dSd/logDp [um2/cm3]"""
        return self.data[self.bin_labels].mul(self.s_multiplier)

    @property
    def dv(self):
        """"""
        return self.dvdlogdp.mul(self.dlogdp)

    @property
    def dvdlogdp(self):
        """Return the histogram in units of dSd/logDp [um3/cm3]"""
        return self.data[self.bin_labels].mul(self.v_multiplier)

    def copy(self):
        """Return a copy of the GenericParticleSizer"""
        return copy.deepcopy(self)

    def dump(self, filepath, ftype='obj'):
        """Save a copy to disk at `filepath`

        :param filepath:
        :param ftype: ['obj', 'json']
        """
        return joblib.dump(self, filepath)

    @property
    def scan_stats(self):
        """"""
        return self.data[list(set(self.data.columns) - set(self.bin_labels))]

    def _subselect_frame(self, df, dmin=0., dmax=1e3):
        """Sub-select the dataframe and return the dataframe

        Units are in microns.
        """
        # generate an array of factors to multiply by
        factors = np.zeros(self.bins.shape[0])

        for i, b in enumerate(self.bins):
            if b[0] >= dmin and b[-1] <= dmax:
                factors[i] = 1
            if dmin >= b[0] and dmin < b[-1]:
                factors[i] = (b[-1] - dmin) / (b[-1] - b[0])
            if dmax >= b[0] and dmax < b[-1]:
                factors[i] = (dmax - b[0]) / (b[-1] - b[0])

        # copy the dataframe
        cpy = df.copy()

        # multiply the histogram by the factors we just calculated
        cpy[self.bin_labels] = cpy[self.bin_labels].mul(factors)

        return cpy

    def stats(self, weight='number'):
        """
        """
        return

    def integrate(self, weight='number', dmin=0., dmax=1., **kwargs):
        """
        """
        rho = kwargs.pop("rho", 1.65)

        weight = weight.lower()

        assert(weight in ['number', 'surface', 'volume', 'mass']), "Invalid `weight`"

        if weight == 'number':
            df = self.dn
        elif weight == 'surface':
            df = self.ds
        elif weight == 'volume':
            df = self.dv
        else:
            df = self.dv * rho

        # subsample the data
        df = self._subselect_frame(df, dmin=dmin, dmax=dmax)

        return df[self.bin_labels].sum(axis=1)

    def slice(self, start=None, end=None, inplace=False):
        """"""

        return

    def resample(self, rs, inplace=False):
        """"""
        return


class SMPS(GenericParticleSizer):
    def __init__(self, **kwargs):
        super(SMPS, self).__init__(**kwargs)

    def __repr__(self):
        return "<SMPS>"


class AlphasenseOpcN2(GenericParticleSizer):
    """
    The Alphasense OPC-N2 is a consumer-grade optical particle counter
    that uses a 658 nm red laser to count and size individual particles. It
    can count particles between 380-17,500 nm.
    """
    def __init__(self, **kwargs):

        # set the bins to be the default bins that
        # Alphasense gives in their datasheet/spreadsheet
        # units are in microns
        bins = kwargs.pop("bins", np.array([
            [0.38, 0.46, 0.54,],
            [0.54, 0.66, 0.78],
            [0.78, 0.915, 1.05],
        ]))

        # default units should be in microns
        dp_units = kwargs.pop("dp_units", "um")

        super(AlphasenseOpcN2, self).__init__(bins=bins, dp_units=dp_units, **kwargs)

    @property
    def bin_weights(self):
        """Allow this to have a setter and a default of 1 for each bin
        """
        return

    def mass(self, dmax=1., rho=1.65):
        """Calculate the particle mass loading between the minimum size cutoff
        and dmax (units are in microns).

        The calculation is done using a discrete integration of each bin.

        :param dmax: particle diameter in units of microns
        :param rho: density in units of g/cm3
        """
        return


class POPS(GenericParticleSizer):
    """
    The Portable Optical Particle Spectrometer is based on the design of R. Gao
    and is currently manufactured by Handix Scientific. It is a research-grade
    OPC that can see particles well below 200 nm.
    """
    def __init__(self, **kwargs):
        bins = kwargs.pop('bins', np.array([]))

        # default units should be microns
        dp_units = kwargs.pop("dp_units", "um")

        super(POPS, self).__init__(bins=bins, dp_units=dp_units, **kwargs)
