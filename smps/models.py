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


__all__ = [
    "GenericParticleSizer", 
    "SMPS", 
    "AlphasenseOPCN2", 
    "AlphasenseOPCN3", 
    "POPS", 
    "ParticlesPlus", 
    "Grimm11D"
]


class GenericParticleSizer(object):
    """
        GenericParticleSizer is the base class for all
        available particle sizing instruments. It contains
        all of the basic functionality and methods used for
        making calculations.

        It can receive data in several formats including: 
        dN, dNdlogDp.
        
        :param data: Data values indexed by timestamp with 
            bins as columns.
        :type data: DataFrame
        :param bins: A list, ordered respective to the bins, 
            where for each bin a min, midpoint, 
            and max particle size is given.
        :type bins: list of lists
        :param dp_units: Units in which the diameter 
            is measured, defaults to 'um'.
        :type dp_units: {'um','nm'}
    """
    def __init__(self, data, bins, **kwargs):
        self.data = data.copy(deep=True)
        self.bins = bins
        self.meta = kwargs.pop('meta', dict())
        self.bin_labels = kwargs.pop('bin_labels', None)

        self.meta['weight'] = kwargs.pop("weight", "number")
        self.meta['units'] = kwargs.pop("units", "dw/dlogdp")

        self.bin_weights = kwargs.pop("bin_weights", 
                                      np.ones(len(self.bins)))

        # if bin_labels were not provided, try to generate them
        if self.bin_labels is None:
            self.bin_labels = [c for c in self.data.columns if 
                               kwargs.pop("bin_prefix", "bin") in c]

        # if no bin labels exist, raise an error
        assert(len(self.bin_labels) > 0), ("No bin labels have been "
        + "found or specified.")

        # make sure bins is has three columns
        assert(self.bins.shape[1] == 3), "`bins` must be an Nx3 array."

        # we always work with our units in um, so convert
        if kwargs.pop('dp_units', 'um') == 'nm':
            self.bins = self.bins * 1e-3

        # if the data is in a format other than 
        # [weight='number', units='dw/dlogdp']
        if kwargs.pop('fmt', 'dn') == 'dn':
            self.data[self.bin_labels] = (
                self
                .data[self.bin_labels]
                .div(self.dlogdp)
            )

        # multiply everything by the bin weights
        # if ones, this will have no effect
        self.data[self.bin_labels] = (
            self
            .data[self.bin_labels]
            .mul(self.bin_weights)
        )

    @property
    def s_multiplier(self):
        """
        A list of the surface areas corresponding
        to the midpoitns of the respective bins.
        """
        return 4 * np.pi * (self.bins[:, 1]/2)**2

    @property
    def v_multiplier(self):
        """
        A list of the volumes corresponding to the 
        midpoitns of the respective bins.
        """
        return (4./3)*np.pi*(self.bins[:, 1]/2)**3

    @property
    def midpoints(self):
        """
        A list of the midpoints of the respective bins.
        """
        return self.bins[:, 1]

    @property
    def dlogdp(self):
        """
        A list of the differences between the upper 
        and lower bound of each bin on a Log scale.
        """
        return np.log10(self.bins[:, -1]) - np.log10(self.bins[:, 0])

    @property
    def dn(self):
        """
        A DataFrame with the number of particles in 
        each bin at each time.
        """
        return self.dndlogdp.mul(self.dlogdp)

    @property
    def dndlogdp(self):
        """
        The histogram in units of dNd/logDp [#/cm3].
        """
        return self.data[self.bin_labels]

#     @property
#     def dd(self):
#         return self.dddlogdp.mul(self.dlogdp)

    @property
    def dddlogdp(self):
        """
        The histogram in units of dDpd/logDp [um/cm3].
        """
        return self.dndlogdp.mul(self.midpoints)

    @property
    def ds(self):
        """
        A DataFrame with the total surface area of 
        all the particles in each bin at each time.
        """
        return self.dn.mul(self.s_multiplier)

    @property
    def dsdlogdp(self):
        """
        The histogram in units of dSd/logDp [um2/cm3].
        """
        return self.dndlogdp.mul(self.s_multiplier)

    @property
    def dv(self):
        """
        A DataFrame with the total volume of all 
        the particels in each bin at each time.
        """
        return self.dvdlogdp.mul(self.dlogdp)

    @property
    def dvdlogdp(self):
        """
        The histogram in units of dSd/logDp [um3/cm3].
        """
        return self.dndlogdp.mul(self.v_multiplier)

    def copy(self):
        """
        Return a copy of the GenericParticleSizer.
        """
        return copy.deepcopy(self)

    def dump(self, filepath, ):
        """
        Save a copy to disk at `filepath`.

        :param filepath: The path at which to save the file.
        :type filepath: string
        :return: The list of file names in which the data is stored.
        :rtype: list of strings
        """
        return joblib.dump(self, filepath)

    @property
    def scan_stats(self):
        """
        A DataFrame with all the data except for the bin data.
        """
        return self.data[list(set(self.data.columns)
                              - set(self.bin_labels))]

    def _subselect_frame(self, df, dmin=0., dmax=1e3):
        """
        Sub-select the bins that are used to make calculations by
        zeroing out bins which are outside of the range. This
        is similar to a slice, except in the bin-dimension 
        rather than the time-dimension. Edge cases are 
        appropriately weighted.

        Units are in microns.
        
        :param df: A DataFrame with bins as columns
        :type df: DataFrame
        :param dmin: The minimum particle diameter in microns to 
            be included.
        :type dmin: float
        :param dmax: The maximum particle diameter in microns to 
            be included.
        :type dmax: float
        :return: A copy of `df` with the undesired columns
            zeroed out.
        :rtype: DataFrame
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

    def stats(self, weight='number', dmin=0., 
              dmax=1e3, rho=1.65, **kwargs):
        """
        Calculate and return the total number of particles, 
        total surface area, total volume, and total mass 
        between dmin and dmax. In addition, the arithmetic 
        mean (AM), geometric mean (GM), mode (Mode), and geometric 
        standard deviation (GSD, CMD) are calculated and are
        weighted by the input parameter `weight`.
        
        :param weight: The dimension to consider in computing 
            statistics.
        :type weight: {"number", "surface", "volume"}
        :param dmin: The minimum particle diameter in microns, 
            defaults to `0`.
        :type dmin: float
        :param dmax: The maximum particle diameter in microns, 
            defaults to `1e3`.
        :type dmax: float
        :param rho: The particle density in units of g/cc3, 
            defaults to `1.65`.
        :type rho: float
        :return: The statistics for each timestamp.
        :rtype: DataFrame
        """
        # remove the weight from the kwargs
        assert(weight in ["number", "surface", "volume", "mass"])

        # initialize an empty dataframe to hold the results
        res = pd.DataFrame(index=self.data.index)

        # subselect the dataframe to only include 
        # diameters of interest
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
            res.loc[idx, "AM"] = (
                1e3*cpy.dn.mul(self.midpoints).sum(axis=1) 
                / res.loc[idx, "number"]
            )
            res.loc[idx, "GM"] = 1e3*np.exp(
                cpy.dn.mul(np.log(self.midpoints), axis=1).sum(axis=1) 
                / res.loc[idx, "number"]
            )
            rename_cols_dict = dict(zip(
                cpy.dndlogdp.columns, 
                cpy.midpoints
            ))
            res.loc[idx, "Mode"] = 1e3 * (
                cpy
                .dndlogdp
                .rename(columns=rename_cols_dict)
                .idxmax(axis=1)
            )

            tmp = cpy.dn.assign(GM=res.loc[idx, 'GM'].values)
        elif weight == "surface":
            res.loc[idx, "AM"] = 1e3 * (
                cpy
                .ds
                .mul(self.midpoints)
                .sum(axis=1) 
                / res.loc[idx, "surface_area"]
            )
            res.loc[idx, "GM"] = 1e3 * np.exp(
                cpy
                .ds
                .mul(np.log(self.midpoints), axis=1)
                .sum(axis=1) 
                / res.loc[idx, "surface_area"]
            )
            rename_cols_dict = dict(zip(
                cpy.dsdlogdp.columns, 
                cpy.midpoints
            ))
            res.loc[idx, "Mode"] = 1e3 * (
                cpy
                .dsdlogdp
                .rename(columns=rename_cols_dict)
                .idxmax(axis=1)
            )

            tmp = cpy.ds.assign(GM=res.loc[idx, 'GM'].values)
        else:
            res.loc[idx, "AM"] = 1e3 * (
                cpy
                .dv
                .mul(self.midpoints)
                .sum(axis=1) 
                / res.loc[idx, "volume"]
            )
            res.loc[idx, "GM"] = 1e3 * np.exp(
                cpy
                .dv
                .mul(np.log(self.midpoints), axis=1)
                .sum(axis=1) 
                / res.loc[idx, "volume"]
            )
            rename_cols_dict = dict(zip(
                cpy.dvdlogdp.columns, 
                cpy.midpoints))
            res.loc[idx, "Mode"] = 1e3 * (
                cpy
                .dvdlogdp
                .rename(columns=rename_cols_dict)
                .idxmax(axis=1)
            )

            tmp = cpy.dv.assign(GM=res.loc[idx, 'GM'].values)

        # calculate the GSD
        res.loc[idx, "GSD"] = self._gsd(tmp)

        # delete the cpy to free up memory
        del cpy, tmp

        return res

    def integrate(self, weight='number', dmin=0., dmax=1., **kwargs):
        """
        For each timestamp give the total number 
        (or surface, volume, or mass) of all the particles 
        (with a diameter between `dmin` and `dmax`) at that time.
        
        :param weight: The variable to sum over, 
            defaults to 'number'.
        :type weight: {'number', 'surface', 'volume', 'mass'}
        :param dmin: The minimum particle diameter in microns, 
            defaults to `0`.
        :type dmin: float
        :param dmax: The maximum particle diameter in microns, 
            defaults to `1e3`.
        :type dmax: float
        :param rho: The particle density in units of g/cc3, 
            defaults to `1.65`.
        :type rho: float
        :return: The desired values indexed by timestamp.
        :rtype: Pandas Series
        """
        rho = kwargs.pop("rho", 1.65)

        weight = weight.lower()

        assert(weight in ['number', 'surface', 'volume']), \
        "Invalid `weight`"

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
        """
        Slice the data between the start and end dates, inclusive.
        
        :param start: The first timestamp to be contained 
            in the result.
        :type start: string, Pandas Timestamp, or numpy DateTime64
        :param end: The last timestamp to be contained in the 
            result.
        :type end: string, Pandas Timestamp, or numpy DateTime64
        :param inplace: If False a new object is returned, 
            defaults to False.
        :type inplace: bool
        :return: The data restricted to the given time interval.
        :rtype: an object the same type as self
        """
        if inplace:
            self.data = self.data[start:end]
        else:
            cpy = self.copy()

            cpy.data = cpy.data[start:end]

            return cpy
        return

    def resample(self, rs, inplace=False):
        """
        Resample on the timeseries.

        :param rs: The resample period using normal 
            Pandas conventions.
        :type rs: string
        :param inplace: If True, update the current instance, 
            otherwise create a new one and return it.
        :type inpace: bool
        :return: A timeseries that has been resampled with 
            period `rs`.
        :rtype: an object from the same class as `self`

        Example::
            m.resample("15min", True)

        """
        obj_cols = (
            self
            .data
            .select_dtypes(include=['object'])
            .resample(rs)
            .first()
        )
        num_cols = (
            self
            .data
            .select_dtypes(exclude=['object'])
            .resample(rs)
            .mean()
        )

        # re-merge the two dataframes
        merged = pd.merge(num_cols, obj_cols, left_index=True, 
                          right_index=True, how='outer')

        if inplace:
            self.data = merged
        else:
            cpy = self.copy()
            cpy.data = merged

            return cpy

        return

    def _gsd(self, df):
        """
        Private function to calculate the geometric 
        standard deviation for a lognormal distribution.
        
        :param df: A DataFrame with a `GM` column holding 
            the geometric mean and columns of bins. 
            In each row the bin value will act as the 
            weight on the midpoint value for that bin, 
            found in `self.midpoints`. 
        :type df: DataFrame
        :return: A column holding the geometric mean for each row.
        :rtype: Pandas Series
        """

        # calculate the GM in units of microns
        gm = df['GM'] * 1e-3

        # drop GM. At this point, the row will have the 
        # histogram and nothing else
        df = df.drop("GM", axis=1)

        # calculate the geometric standard deviation
        gsd = np.exp(
            np.sqrt(
                df
                 .mul(
                     (df*0)
                     .add(np.log(self.midpoints), axis=1)
                     .sub(np.log(gm), axis=0)
                     **2
                 )
                 .sum(axis=1)
                / (df.sum(axis=1))
            )
        )

        return gsd

    
class SMPS(GenericParticleSizer):
    def __init__(self, **kwargs):
        super(SMPS, self).__init__(fmt='dndlogdp', dp_units='nm', 
                                   **kwargs)
    def __repr__(self):
        return "<SMPS>"


class AlphasenseOPCN2(GenericParticleSizer):
    """
    The Alphasense OPC-N2 is a consumer-grade optical 
    particle counter that uses a 658 nm red laser to count 
    and size individual particles. It can count particles 
    between 380-17,500 nm. Read more `here
    <https://www.alphasense.com/products/
    optical-particle-counter/>`__.
    """
    def __init__(self, **kwargs):
        bb = np.array([0.38, 0.54, 0.78, 1.05, 1.34, 1.59, 2.07, 3., 
                       4., 5., 6.5, 8., 10., 12., 14., 16., 17.5])

        bins = kwargs.pop("bins", make_bins(boundaries=bb))
        fmt = kwargs.pop("fmt", "dn")
        bin_labels = kwargs.pop(
            "bin_labels", 
            ["bin{}".format(i) for i in range(bins.shape[0])]
        )

        super(AlphasenseOPCN2, self).__init__(
            bins=bins, 
            fmt=fmt, 
            bin_labels=bin_labels, 
            **kwargs
        )


class AlphasenseOPCN3(GenericParticleSizer):
    """
    The Alphasense OPC-N3 is a consumer-grade optical
    particle counter that uses a 658 nm red laser to count 
    and size individual particles. It can count particles 
    between 380-40,000 nm. Read more `here <https://www.
    alphasense.com/products/optical-particle-counter/>`__.
    """
    def __init__(self, **kwargs):
        bb = np.array([0.35, 0.46, 0.66, 1., 1.3, 1.7, 2.3, 3.,
                            4., 5.2, 6.5, 8., 10., 12., 14., 16., 18.,
                            20., 22., 25., 28., 31., 34., 37., 40.])

        bins = kwargs.pop("bins", make_bins(boundaries=bb))
        fmt = kwargs.pop("fmt", "dn")
        bin_labels = kwargs.pop(
            "bin_labels", 
            ["bin{}".format(i) for i in range(bins.shape[0])]
        )

        super(AlphasenseOPCN3, self).__init__(
            bins=bins, 
            fmt=fmt, 
            bin_labels=bin_labels, 
            **kwargs
        )


class POPS(GenericParticleSizer):
    """
    The Portable Optical Particle Spectrometer is based on 
    the design of R. Gao and is currently manufactured by 
    Handix Scientific. It is a research-grade OPC that can see 
    particles as small as 132 nm. Read more `here 
    <http://www.handixscientific.com/pops>`__.
    """
    def __init__(self, **kwargs):
        # default left and right bounds according to Handix
        bb = 1e-3*np.array([132, 144, 158, 174, 191, 209, 229, 270, 324, 
                            473, 594, 1009, 1294, 1637, 2148, 2864, 3648])

        # make the 3xn array of bins using the utility 
        # function make_bins()
        bins = kwargs.pop('bins', make_bins(boundaries=bb))

        #
        fmt = kwargs.pop("fmt", "dn")

        super(POPS, self).__init__(
            bins=bins, 
            fmt=fmt, 
            bin_labels=bin_labels, 
            **kwargs
        )


class ParticlesPlus(GenericParticleSizer):
    """
    The Particles Plus is a consumer-grade optical particle counter
    that uses a red laser to count and size individual particles. It
    can count particles between 300-10,000 nm. Read more `here 
    <https://particlesplus.com/ambient-air-monitoring/>`__.
    """
    def __init__(self, **kwargs):
        bb = np.array([0.3, 0.5, 0.7, 1., 2.5, 10.0, 10.1])

        bins = kwargs.pop("bins", make_bins(boundaries=bb))
        fmt = kwargs.pop("fmt", "dn")
        bin_labels = kwargs.pop(
            "bin_labels", 
            ["bin{}".format(i) for i in range(bins.shape[0])]
        )

        super(ParticlesPlus, self).__init__(
            bins=bins, 
            fmt=fmt, 
            bin_labels=bin_labels, 
            **kwargs
        )


class Grimm11D(GenericParticleSizer):
    """
    The Grimm11D is an optical particle counter that uses a blue 
    laser to count and size individual particles. It can count 
    particles between 250-35,000 nm. Read more `here 
    <https://www.grimm-aerosol.com/products-en/dust-monitors/the-
    dust-decoder/11-d/>`__.
    """
    def __init__(self, **kwargs):
        bb = np.array([
            0.25, 0.28, 0.3, 0.35, 0.4, 0.45, 0.5, 0.58, 0.65, 0.7, 
            0.8, 1.0, 1.3, 1.6, 2.0, 2.5, 3.0, 3.5, 4., 5., 6.5, 7.5, 
            8.5, 10., 12.5, 15., 17.5, 20., 25., 30., 32., 35.
        ])

        bins = kwargs.pop("bins", make_bins(boundaries=bb))
        fmt = kwargs.pop("fmt", "dn")
        bin_labels = kwargs.pop(
            "bin_labels", 
            ["bin{}".format(i) for i in range(bins.shape[0])]
        )

        super(Grimm11D, self).__init__(
            bins=bins, 
            fmt=fmt, 
            bin_labels=bin_labels, 
            **kwargs
        )
