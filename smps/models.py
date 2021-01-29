#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import pandas as pd
import math
import copy
import joblib
import json
from .utils import make_bins

__all__ = ["GenericParticleSizer", "SMPS", "AlphasenseOpcN2", "AlphasenseOpcN3", 
            "POPS", "ParticlesPlus", "Grimm11D"]


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

        self.bin_weights = kwargs.pop("bin_weights", np.ones(len(self.bins)))

        # if bin_labels were not provided, try to generate them
        if self.bin_labels is None:
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

        # multiply everything by the bin weights
        # if ones, this will have no effect
        self.data[self.bin_labels] = self.data[self.bin_labels].mul(self.bin_weights)

    @property
    def s_multiplier(self):
        return 4 * np.pi * (self.bins[:, 1]/2)**2

    @property
    def v_multiplier(self):
        return (4./3)*np.pi*(self.bins[:, 1]/2)**3

    @property
    def midpoints(self):
        return self.bins[:, 1]

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
        return self.dndlogdp.mul(self.midpoints)

    @property
    def ds(self):
        """"""
        return self.dn.mul(self.s_multiplier)

    @property
    def dsdlogdp(self):
        """Return the histogram in units of dSd/logDp [um2/cm3]"""
        return self.dndlogdp.mul(self.s_multiplier)

    @property
    def dv(self):
        """"""
        return self.dvdlogdp.mul(self.dlogdp)

    @property
    def dvdlogdp(self):
        """Return the histogram in units of dSd/logDp [um3/cm3]"""
        return self.dndlogdp.mul(self.v_multiplier)

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
        """Sub-select the bins that are used to make calculations. This
        is similar to a slice, except in the bin-dimension rather than the
        time-dimension.

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

    def stats(self, weight='number', dmin=0., dmax=1e3, rho=1.65, **kwargs):
        """Calculate and return the total number of particles, total
        surface area, total volume, and total mass between dmin and dmax. In
        addition, the arithmetic mean (AM), geometric mean (GM), mode (Mode),
        and geometric standard deviation (GSD, CMD) are calculated and are
        weighted by the input parameter `weight`.

        :param dmin: minimum particle diameter in microns
        :param dmax: maximum particle diameter in microns
        :param rho: particle density in units of g/cc3
        """
        # remove the weight from the kwargs
        assert(weight in ["number", "surface", "volume", "mass"])

        # initialize an empty dataframe to hold the results
        res = pd.DataFrame(index=self.data.index)

        # subselect the dataframe to only include diameters of interest
        cpy = self.copy()
        cpy.data = cpy._subselect_frame(cpy.data, **kwargs)

        # drop all rows that are completely NaN
        cpy.data = cpy.data.dropna(how='all')

        # make a shortcut for the index to inject
        idx = cpy.data.index

        # calculate the total number of particles
        res.loc[idx, "number"] = cpy.dn.sum(axis=1)
        res.loc[idx, "surface_area"] = cpy.ds.sum(axis=1)
        res.loc[idx, "volume"] = cpy.dv.sum(axis=1)
        res.loc[idx, "mass"] = cpy.dv.mul(rho).sum(axis=1)

        if weight == "number":
            res.loc[idx, "AM"] = 1e3*cpy.dn.mul(self.midpoints).sum(axis=1) / res.loc[idx, "number"]
            res.loc[idx, "GM"] = 1e3*np.exp(cpy.dn.mul(np.log(self.midpoints), axis=1).sum(axis=1) / res.loc[idx, "number"])
            res.loc[idx, "Mode"] = cpy.dndlogdp.apply(lambda x: 1e3*cpy.midpoints[cpy.dndlogdp.columns.get_loc(x.idxmax())], axis=1)

            tmp = cpy.dn.assign(GM=res.loc[idx, 'GM'].values)
        elif weight == "surface":
            res.loc[idx, "AM"] = 1e3 * cpy.ds.mul(self.midpoints).sum(axis=1) / res.loc[idx, "surface_area"]
            res.loc[idx, "GM"] = 1e3 * np.exp(cpy.ds.mul(np.log(self.midpoints), axis=1).sum(axis=1) / res.loc[idx, "surface_area"])
            res.loc[idx, "Mode"] = cpy.dsdlogdp.apply(lambda x: 1e3*cpy.midpoints[cpy.dsdlogdp.columns.get_loc(x.idxmax())], axis=1)

            tmp = cpy.ds.assign(GM=res.loc[idx, 'GM'].values)
        else:
            res.loc[idx, "AM"] = 1e3 * cpy.dv.mul(self.midpoints).sum(axis=1) / res.loc[idx, "volume"]
            res.loc[idx, "GM"] = 1e3 * np.exp(cpy.dv.mul(np.log(self.midpoints), axis=1).sum(axis=1) / res.loc[idx, "volume"])
            res.loc[idx, "Mode"] = cpy.dvdlogdp.apply(lambda x: 1e3*cpy.midpoints[cpy.dvdlogdp.columns.get_loc(x.idxmax())], axis=1)

            tmp = cpy.dv.assign(GM=res.loc[idx, 'GM'].values)

        # calculate the GSD
        res.loc[idx, "GSD"] = tmp.apply(self._gsd, axis=1)

        # delete the cpy to free up memory
        del cpy, tmp

        return res

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
        """Slice the data between the start and end dates
        """
        if inplace:
            self.data = self.data[start:end]
        else:
            cpy = self.copy()

            cpy.data = cpy.data[start:end]

            return cpy
        return

    def resample(self, rs, inplace=False):
        """Resample on the timeseries.

        :param rs: resample period using normal Pandas conventions
        :param inplace: if True, update the current instance, else create a new one and return

        Example:
        >>> m.resample("15min", True)

        """
        obj_cols = self.data.select_dtypes(include=['object']).resample(rs).first()
        num_cols = self.data.select_dtypes(exclude=['object']).resample(rs).mean()

        # re-merge the two dataframes
        merged = pd.merge(num_cols, obj_cols, left_index=True, right_index=True, how='outer')

        if inplace:
            self.data = merged
        else:
            cpy = self.copy()
            cpy.data = merged

            return cpy

        return

    def _gsd(self, row):
        """Private function to calculate the geometric standard deviation for a
        lognormal distribution. The equation used is:

        log...
        """
        # find the row that the geometric mean is in
        idx = row.index.isin(["GM"])

        # calculate the GM in units of microns
        gm = row.loc[idx]['GM'] * 1e-3

        # grab the row with the exception of the GM
        # at this point, the row will have the histogram and nothing else
        row = row.loc[~idx]

        # calculate the geometric standard deviation
        gsd = np.exp(np.sqrt((row.mul((np.log(self.midpoints) - np.log(gm))**2).sum()) / (row.sum() - 1)))

        return gsd


class SMPS(GenericParticleSizer):
    def __init__(self, **kwargs):
        super(SMPS, self).__init__(fmt='dndlogdp', dp_units='nm', **kwargs)

    def __repr__(self):
        return "<SMPS>"


class AlphasenseOpcN2(GenericParticleSizer):
    """
    The Alphasense OPC-N2 is a consumer-grade optical particle counter
    that uses a 658 nm red laser to count and size individual particles. It
    can count particles between 380-17,500 nm.
    """
    def __init__(self, **kwargs):
        bb = np.array([0.38, 0.54, 0.78, 1.05, 1.34, 1.59, 2.07, 3.,
                            4., 5., 6.5, 8., 10., 12., 14., 16., 17.5])

        bins = kwargs.pop("bins", make_bins(boundaries=bb))

        super(AlphasenseOpcN2, self).__init__(bins=bins, fmt='dn', **kwargs)


class AlphasenseOpcN3(GenericParticleSizer):
    """
    The Alphasense OPC-N3 is a consumer-grade optical particle counter
    that uses a 658 nm red laser to count and size individual particles. It
    can count particles between 380-40,000 nm.
    """
    def __init__(self, **kwargs):
        bb = np.array([0.35, 0.46, 0.66, 1., 1.3, 1.7, 2.3, 3.,
                            4., 5.2, 6.5, 8., 10., 12., 14., 16., 18.,
                            20., 22., 25., 28., 31., 34., 37., 40.])

        bins = kwargs.pop("bins", make_bins(boundaries=bb))

        super(AlphasenseOpcN3, self).__init__(bins=bins, fmt='dn', **kwargs)


class POPS(GenericParticleSizer):
    """
    The Portable Optical Particle Spectrometer is based on the design of R. Gao
    and is currently manufactured by Handix Scientific. It is a research-grade
    OPC that can see particles as small as 132 nm.
    """
    def __init__(self, **kwargs):
        # default left and right bounds according to Handix
        bb = 1e-3*np.array([132, 144, 158, 174, 191, 209, 229, 270, 324, 473, 594,
                        1009, 1294, 1637, 2148, 2864, 3648])

        # make the 3xn array of bins using the utility function make_bins()
        bins = kwargs.pop('bins', make_bins(boundaries=bb))

        super(POPS, self).__init__(bins=bins, **kwargs)


class ParticlesPlus(GenericParticleSizer):
    """
    The Particles Plus is a consumer-grade optical particle counter
    that uses a red laser to count and size individual particles. It
    can count particles between 300-10,000 nm.
    """
    def __init__(self, **kwargs):
        bb = np.array([0.3, 0.5, 0.7, 1., 2.5, 10.0, 10.1])

        bins = kwargs.pop("bins", make_bins(boundaries=bb))

        super(ParticlesPlus, self).__init__(bins=bins, fmt='dn', **kwargs)


class Grimm11D(GenericParticleSizer):
    """
    The Grimm11D is an optical particle counter
    that uses a blue laser to count and size individual particles. It
    can count particles between 250-35,000 nm.
    """
    def __init__(self, **kwargs):
        bb = np.array([0.25, 0.28, 0.3, 0.35, 0.4, 0.45, 0.5, 0.58, 0.65, 0.7, 0.8,
                        1.0, 1.3, 1.6, 2.0, 2.5, 3.0, 3.5, 4., 5., 6.5, 7.5, 8.5, 10.,
                        12.5, 15., 17.5, 20., 25., 30., 32., 35.])

        bins = kwargs.pop("bins", make_bins(boundaries=bb))

        super(Grimm11D, self).__init__(bins=bins, fmt='dn', **kwargs)