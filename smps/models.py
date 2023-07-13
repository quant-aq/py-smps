#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import pandas as pd
import math
import copy
import sys
import joblib
import json
import logging
from .utils import make_bins

logging.basicConfig(stream=sys.stdout, level=logging.WARN)

__all__ = [
    "GenericParticleSizer", 
    "SMPS", 
    "AlphasenseOPCN2", 
    "AlphasenseOPCN3", 
    "ModulairPM", 
    "Modulair", 
    "POPS", 
    "ParticlesPlus", 
    "Grimm11D"
]


class ValidationError(Exception):
    def __init__(self, msg):
        self.msg = msg


class GenericParticleSizer(object):
    """
    The base class for a Generic Size-Resolving Particle Instrument.
    
    Parameters
    ----------
    data : pd.DataFrame
            A dataframe containing the binned particle data, timestamp column, 
            and optional metadata columns.
    bins :  array-like
            A 3xn array defining the left, mid, and right boundaries for 
            each particle size bin.
    dp_units : str, deafult='um'
        Define the units of particle bin sizes - one of ('um', 'nm')
    meta : dict, default=None
        An optional dictionary containing meta information about the instrument.
    weight : str, default='number'
        The moment of the model. One of ('number', 'surface', 'volume').
    units : str, default='dw/dlogdp'
        The default data format in AIM-compatible string format.
    bin_weights : array-like, default=ones
        An optional array to set bin weights (e.g., counting efficiency). By 
        default, all bin weights are 1.
    bin_prefix : str, default='bin'
        An optional kwarg to name the bins if not already named.
    fmt : str, default='dn' 
        Optional kwarg describing the default data format. One of ('dn', 'dndlogdp').
    
    See Also
    --------
    smps.models.Modulair
    smps.models.ModulairPM
    smps.models.SMPS
    smps.models.Grimm11D
    smps.models.POPS
    smps.models.ParticlesPlus
    smps.models.AlphasenseOPCN2
    smps.models.AlphasenseOPCN3

    Examples
    --------
    
    Initialize a base object with data and bins:
    
    >>> obj = smps.models.GenericParticleSizer(data=df, bins=bins)
    
    """
    def __init__(self, data, bins, **kwargs):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)

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
        Per particle surface area representing the mean particle in each bin.
        """
        return 4 * np.pi * (self.bins[:, 1]/2)**2

    @property
    def v_multiplier(self):
        """
        Per particle volume representing the mean particle in each bin.
        """
        return (4./3)*np.pi*(self.bins[:, 1]/2)**3

    @property
    def midpoints(self):
        """
        Midpoint particle diameter in each bin.
        """
        return self.bins[:, 1]

    @property
    def dlogdp(self):
        """
        Log difference between the upper and lower bound of each bin.
        """
        return np.log10(self.bins[:, -1]) - np.log10(self.bins[:, 0])

    @property
    def dn(self):
        """
        The number concentration per bin. Can otherwise be represented as :math:`dN/dD_p`. Units 
        are pp/cm3.
        """
        return self.dndlogdp.mul(self.dlogdp)

    @property
    def dndlogdp(self):
        """
        The log normalized number concentration per bin. Can otherwise be represented as 
        :math:`dN/dlogD_p`. Units are pp/cm3.
        """
        return self.data[self.bin_labels]

    @property
    def dddlogdp(self):
        """
        Can otherwise be represented as :math:`dD/dlogD_p`.
        """
        return self.dndlogdp.mul(self.midpoints)

    @property
    def ds(self):
        """
        Surface area per bin.
        """
        return self.dn.mul(self.s_multiplier)

    @property
    def dsdlogdp(self):
        """
        Log normalized surface area per bin. Otherwise represented as :math:`dS/dlogDp`.
        Units are: µm2/cm3.
        """
        return self.dndlogdp.mul(self.s_multiplier)

    @property
    def dv(self):
        """
        Volume per bin.
        """
        return self.dvdlogdp.mul(self.dlogdp)

    @property
    def dvdlogdp(self):
        """
        Log normalized volume per bin. Can otherwise be represented as
        :math:`dV/dlogDp`. Units are µm3/cm3.
        """
        return self.dndlogdp.mul(self.v_multiplier)

    @property
    def scan_stats(self):
        """
        Return the non-binned data.
        """
        return self.data[list(set(self.data.columns)
                              - set(self.bin_labels))]
        
    def copy(self):
        """
        Return a copy of the GenericParticleSizer.
        
        Examples
        --------
        >>> obj = GenericParticleSize(data=df)
        >>> cpy = obj.copy()
        """
        return copy.deepcopy(self)

    def dump(self, filepath):
        """
        Save a copy of the model to disk.
        
        Parameters
        ----------
        filepath : str
            The filepath to save the file at.
        
        Returns
        -------
        filepath : list
            The list of filepaths where data was saved.
        
        Examples
        --------
        
        >>> model.dump(filepath="path-to-file.sav")
        
        """
        return joblib.dump(self, filepath)

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
    
    def _subselect_bins(self, **kwargs):
        """Return an array of multipliers corresponding to the 
        percentage of a given bin to use in a calculation based 
        on the size cutoffs defined.
        
        Args:
            bins (array): Array of bins
            dmin (float): Minimum particle diameter
            dmax (float): Maximum particle diameter
        """
        bins = kwargs.pop("bins")
        dmin = kwargs.pop("dmin")
        dmax = kwargs.pop("dmax")
        
        # Create a placeholder for the factors
        factors = np.zeros(bins.shape[0])
        
        for i, b in enumerate(bins):
            if b[0] >= dmin and b[-1] <= dmax:
                factors[i] = 1.
            
            if dmin >= b[0] and dmin < b[-1]:
                factors[i] = (b[-1] - dmin) / (b[-1] - b[0])
            
            if dmax >= b[0] and dmax < b[-1]:
                factors[i] = (dmax - b[0]) / (b[-1] - b[0])
        
        return factors

    def stats(self, weight='number', dmin=0., dmax=1e3, rho=1.65, **kwargs):
        """
        Compute and return the aerosol size distribution statistics.
        
        Calculate the total number of particles, surface area, 
        volume, mass, arithmetic mean particle diameter, 
        geometric mean particle diameter, and geometric 
        standard deviation between any two particle diameters.
        
        Parameters
        ----------
        weight: str, default='number'
            The moment to use to compute the statistics. Should be one of 
            ('number', 'surface', 'volume', 'mass').
        dmin : float, default=0.0
            The minimum particle diameter to consider. Units are in µm.
        dmax : float, default=1000.0
            The maximum particle diameter to consider. Units are in µm.
        rho : float, default=1.65
            The particle density in units of g/cm3.
        
        Returns
        -------
        data : pd.DataFrame
            A dataframe containing the computed statistics at each point in time.
        
        Examples
        --------
        
        Compute the number-weighted stats
        
        >>> obj.stats(weight='number')
        
        Compute the volume-weighted stats for PM2.5
        
        >>> obj.stats(weight='volume', dmin=0, dmax=2.5)
        
        """
        # remove the weight from the kwargs
        assert(weight in ["number", "surface", "volume", "mass"])

        # initialize an empty dataframe to hold the results
        res = pd.DataFrame(index=self.data.index)

        # subselect the dataframe to only include 
        # diameters of interest
        cpy = self.copy()
        cpy.data = cpy._subselect_frame(cpy.data, dmin=dmin, dmax=dmax, **kwargs)

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
        
        if weight == "mass":
            res.loc[idx, "AM"] = res.loc[idx, "AM"] * rho
            res.loc[idx, "GM"] = res.loc[idx, "GM"] * rho
            res.loc[idx, "Mode"] = res.loc[idx, "Mode"] * rho
        
        
        # calculate the GSD
        res.loc[idx, "GSD"] = self._gsd(tmp)

        # delete the cpy to free up memory
        del cpy, tmp

        return res

    def integrate(self, weight='number', dmin=0., dmax=1., **kwargs):
        """
        Integrate by weight for each record in time.
        
        Compute the total number of particles, surface area, volume, or 
        mass between any two particle diameters. Correct for hygroscopic 
        growth by defining a hygroscopic growth factor, kappa.
        
        Parameters
        ----------
        weight : str, default='number'
            The moment/variable to integrate. One of ('number', 'surface', 'volume', 'mass').
        dmin : float, default=0.0
            The minimum particle diameter in µm.
        dmax : float, default=1.
            The maximum particle diameter in µm.
        rho : float or callable
            The density as a float or a callable function which takes particle diameter as its 
            only argument and returns the density at that diameter.
        kappa : float or callable
            The kappa growth factor as a float or a callable function which takes particles diameter 
            as its only argument and returns the kappa value at that diameter.
        rh : str
            If kappa is defined, rh is the column of data corresponding to the relative humidity 
            which is required to correct for hygroscopic growth.
        
        Returns
        -------
        data : pd.Series
            A series containing the integrated values.
        
        Examples
        --------
        
        Compute the total particle concentration under 1 µm:
        
        >>> obj.integrate(weight='number', dmin=0., dmax=1.)
        
        Compute PM1
        
        >>> obj.integrate(weight='mass', dmin=0., dmax=1.)
        
        Compute PM2.5 assuming a kappa=0.3
        
        >>> obj.integrate(weight='mass', dmin=0., dmax=2.5, kappa=0.3, rh="sample_rh")
        
        """
        # Convert the density to a callable if not already
        rho = kwargs.pop("rho", 1.65)
        if not callable(rho):
            _rho = lambda x: rho
        else:
            _rho = rho

        assert(weight in ['number', 'surface', 'volume', 'mass']), \
                "Invalid `weight`"
                
        # Check for kappa
        kappa = kwargs.pop("kappa", None)
        if kappa:
            if not callable(kappa):
                _kappa = lambda x: kappa
            else:
                _kappa = kappa
                
            # Check for a column with rh
            rh_c = kwargs.pop("rh", None)
            if not rh_c:
                raise AttributeError("`rh` is a required column when correcting for hygroscopic growth")
            
            # Check to make sure the rh column exists
            if rh_c not in self.data.columns:
                raise ValidationError(f"No column for {rh_c} was found in `data`")
            
            # Check to see if there are rh values below 1% or above 95%
            count_low_rh = self.data[self.data[rh_c] <= 1.0].shape[0]
            if count_low_rh > 0:
                self.logger.warning(f"There are {count_low_rh} values below 1% RH in your data which may cause problems when correcting for hygroscopic growth.")
            
            count_high_rh = self.data[self.data[rh_c] >= 95.0].shape[0]
            if count_high_rh > 0:
                self.logger.warning(f"There are {count_high_rh} values above 95% RH in your data which may cause problems when correcting for hygroscopic growth.")
        else:
            _kappa = None
            
        # If kappa is not set, compute the integrated values and return
        if not _kappa:
            if weight == 'number':
                df = self.dn
            elif weight == 'surface':
                df = self.ds
            elif weight == 'volume':
                df = self.dv
            else:
                df = self.dv * [_rho(dp) for dp in self.midpoints]
            
            # Subsample the data and return
            df = self._subselect_frame(df, dmin=dmin, dmax=dmax)
            
            return df[self.bin_labels].sum(axis=1)
             
        def compute_integration_by_row(row, weight):
            # Recompute the bins
            _bins = self.bins.copy()
            
            # Get rh
            rh = row[rh_c]
            
            # Compute kappa for each (wet) diameter
            k = [_kappa(dp) for dp in _bins[:, 1]]
            
            for i in range(_bins.shape[0]):
                for j in range(3):
                    _bins[i, j] = _bins[i, j] / (1 + k[i] * (rh / (100. - rh)))**(1./3.)
                    
            # Extract the midpoints
            midpoints = _bins[:, 1]
            
            # Keep only the data for the bins that are relevant for this record
            histogram = row[self.bin_labels] * self._subselect_bins(dmin=dmin, dmax=dmax, bins=_bins)
            
            # Return the histogram depending on which weight we're using
            if weight == 'surface':
                hgram = histogram.mul((4 * np.pi * (_bins[:, 1]/2)**2))
            elif weight == 'volume':
                hgram = histogram.mul((4./3)*np.pi*(_bins[:, 1]/2)**3)
            elif weight == 'mass':
                hgram = histogram.mul((4./3)*np.pi*(_bins[:, 1]/2)**3).mul([_rho(dp) for dp in _bins[:, 1]])
            else:
                hgram = histogram
                
            return hgram.sum()
        
        # Create a dataframe that is dn + humidity so that we are always sending the correct 
        # format to the integration computation 
        cpy = self.dn
        cpy.loc[:, rh_c] = self.data[rh_c]
        
        s = cpy.apply(compute_integration_by_row, axis=1, weight=weight)
        
        return s

    def slice(self, start=None, end=None, inplace=False):
        """
        Slice the data between two datetimes.
        
        Parameters
        ----------
        start : str
            The timestamp (as a string) to begin the slice.
        end : str
            The timestamp (as a string) to end the slice.
        inplace : bool, default=False
            If True, replace the data with the specified slice of the data, 
            otherwise return a copy of the sliced data.
        
        Returns
        -------
        data : pd.DataFrame or None
            If inplace=False, return the dataframe with the sliced data.
        
        Examples
        --------
        >>> obj.slice(start="2023-01-01", end="2023-02-15")
        
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
        Resample the data to a new time base.
        
        Parameters
        ----------
        rs : str
            The resample period using normal pandas convention. Examples include 
            '1h', '1d', '15min', etc.
        inplace : bool, default=False
            If True, modify the data inplace, otherwise return.
        
        Returns
        -------
        data : pd.DataFrame or None
            If inplace=False, return the dataframe with the resampled data.
        
        Examples
        --------
        
        Resample to 1h and modify in place
        
        >>> obj.resample('1h', inplace=True)

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
    """
    The TSI Scanning Mobility Particle Sizer.
    
    The TSI SMPS is a research instrument 
    for counting and sizing submicron aerosol.
    
    """
    def __init__(self, **kwargs):
        super(SMPS, self).__init__(fmt='dndlogdp', dp_units='nm', **kwargs)


class AlphasenseOPCN2(GenericParticleSizer):
    """
    The Alphasense OPC-N2.
    
    The Alphasense OPC-N2 is a consumer-grade OPC
    that uses a 658 nm red laser to count 
    and size particles between 380 - 17,500 nm. Read more `here
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
    The Alphasense OPC-N3.
    
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


class ModulairPM(GenericParticleSizer):
    """
    The QuantAQ MODULAIR-PM.
    
    QuantAQ's MODULAIR-PM sensors use the Alphasense OPC-N3,
    a consumer-grade optical particle counter that uses a 658 nm 
    red laser to count and size individual particles. It can 
    count particles between 380-40,000 nm. Read more 
    `here <https://www.alphasense.com/products/optical-particle-counter/>`__.
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

        super(ModulairPM, self).__init__(
            bins=bins, 
            fmt=fmt, 
            bin_labels=bin_labels, 
            **kwargs
        )


class Modulair(GenericParticleSizer):
    """
    The QuantAQ MODULAIR.
    
    QuantAQ's MODULAIR sensors use the Alphasense OPC-N3,
    a consumer-grade optical particle counter that uses a 658 nm 
    red laser to count and size individual particles. It can 
    count particles between 380-40,000 nm. However, because the last
    few bins (bin18-bin24) collect negligible data, they are removed to conserve 
    memory, which means the MODULAIR data differs slightly from the MODULAIR-PM data.
    Read more about the Alphasense OPC-N3
    `here <https://www.alphasense.com/products/optical-particle-counter/>`__.
    """
    def __init__(self, **kwargs):
        bb = np.array([0.35, 0.46, 0.66, 1., 1.3, 1.7, 2.3, 3.,
                            4., 5.2, 6.5, 8., 10., 12., 14., 16., 18.,
                            20., 22.,])

        bins = kwargs.pop("bins", make_bins(boundaries=bb))
        fmt = kwargs.pop("fmt", "dn")
        bin_labels = kwargs.pop(
            "bin_labels", 
            ["bin{}".format(i) for i in range(bins.shape[0])]
        )

        super(Modulair, self).__init__(
            bins=bins, 
            fmt=fmt, 
            bin_labels=bin_labels, 
            **kwargs
        )


class POPS(GenericParticleSizer):
    """
    The Handix Scientific Portable Optical Particle Spectrometer (POPS).
    
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
    The Particles Plus OPC.
    
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
    The GRIMM 11D.
    
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
