#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import pandas as pd
import numpy as np
import math
import copy
import joblib

from .utils import _get_bin_count, _get_linecount, RENAMED_COLUMNS
from .plots import heatmap


class SMPS(object):
    """
        Build an SMPS object from the raw data and bin information. It is assumed
        that data fed into the object is in dN/dlogDp format and any other format
        will lead to errors.
    """
    def __init__(self, data, bins, bin_labels, dp_units='nm', meta=None):
        """
        
        dp_units : ['nm', 'um']
        """
        self.raw = data
        self.meta = meta
        self.bins = bins

        # Convert all diameters to microns
        if dp_units == 'nm':
            self.bins /= 1000.

        self.midpoints = bins[:, 1]
        self.bin_labels = bin_labels

        self.s_multiplier = self.midpoints**2 * np.pi
        self.v_multiplier = self.midpoints**3 * (np.pi/6.)

    def copy(self):
        """Return a copy of the SMPS instance."""
        return copy.deepcopy(self)

    def dump(self, filepath):
        """Save the SMPS object to disk"""
        return joblib.dump(self, filepath)

    @property
    def dlogdp(self):
        """Return the bin midpoints"""
        return np.log10(self.bins[:, -1]) - np.log10(self.bins[:, 0])

    @property
    def dndlogdp(self):
        """Return dN/dlogDp in units of [#/um*cm3]"""
        return self.raw[self.bin_labels]

    @property
    def dddlogdp(self):
        """Return dDp/dlogDp in units of [#/cm3]"""
        return self.dndlogdp.mul(self.midpoints)

    @property
    def dsdlogdp(self):
        """Return dSdlogDp"""
        return self.dndlogdp.mul(self.s_multiplier)

    @property
    def dvdlogdp(self):
        """Return dVdlogDp"""
        return self.dndlogdp.mul(self.v_multiplier)

    @property
    def dn(self):
        """Return dN in units of [#/cm3]"""
        return self.dndlogdp.mul(self.dlogdp)

    @property
    def dd(self):
        """Return some make believe dDp"""
        return self.dddlogdp.mul(self.dlogdp)

    @property
    def ds(self):
        """Return dS un units of [um2/cm3]"""
        return self.dsdlogdp.mul(self.dlogdp)

    @property
    def dv(self):
        """Returns dV in units of [um3/cm3]"""
        return self.dvdlogdp.mul(self.dlogdp)

    def stats(self, weight='number', rho=1.65):
        """Stats can be weighted by: [number, surface_area, volume].

        Columns:
            Total Number: [cm-3]
            Total Surface Area: [um2 cm-3]
            Total Volume: [um3 cm-3]
            Total Mass: [ug m-3]
            Mean: [nm]

        Stats returned include: Total, GM, GSD, Mean, Median, CMD"""
        res = pd.DataFrame()

        res['Total Number'] = self.dn.sum(axis=1)
        res['Total Surface Area'] = self.ds.sum(axis=1)
        res['Total Volume'] = self.dv.sum(axis=1)
        res['Total Mass'] = self.dv.sum(axis=1)*rho

        if weight == 'number':
            res['Mean'] = 1e3 * self.dn.mul(self.midpoints).sum(axis=1) / res['Total Number']
            res['GM'] = 1e3 * np.exp(self.dn.mul(np.log(self.midpoints), axis=1).sum(axis=1) / res['Total Number'])

            # Create a tmp column with just inner part of the GSD calculation
            tmp = self.dn.assign(GM=res['GM'].values)

            res['GSD'] = tmp.apply(self._gsd, axis=1)
        elif weight == 'surface_area': # 1e3 is to convert to um from nm
            res['Mean'] = 1e3 * self.ds.mul(self.midpoints).sum(axis=1) / res['Total Surface Area']
            res['GM'] = 1e3 * np.exp(self.ds.mul(np.log(self.midpoints), axis=1).sum(axis=1) / res['Total Surface Area'])

            # Create a tmp column with just inner part of the GSD calculation
            tmp = self.ds.assign(GM=res['GM'].values)

            res['GSD'] = tmp.apply(self._gsd, axis=1)
        elif weight == 'volume':
            res['Mean'] = 1e3 * self.dv.mul(self.midpoints).sum(axis=1) / res['Total Volume']
            res['GM'] = 1e3 * np.exp(self.dv.mul(np.log(self.midpoints), axis=1).sum(axis=1) / res['Total Volume'])

            # Create a tmp column with just inner part of the GSD calculation
            tmp = self.dv.assign(GM=res['GM'].values)

            res['GSD'] = tmp.apply(self._gsd, axis=1)
        else:
            raise Exception("Invalid parameter for weight.")

        # Clear up some memory
        del tmp

        return res

    def _gsd(self, row):
        """Calculate and return the GSD for a given row.
        """
        gm_idx = row.index.isin(['GM'])
        gm = row.loc[gm_idx]['GM'] * 1e-3
        row = row.loc[~gm_idx]

        return np.exp(np.sqrt((row.mul((np.log(self.midpoints) - np.log(gm))**2).sum())/(row.sum())))

    @property
    def scan_stats(self):
        """Return the scan meta information for each row as a DataFrame."""
        cols_to_get = set(self.raw.columns) - set(self.bin_labels)

        return self.raw[list(cols_to_get)]

    def integrate(self, dmin, dmax, weight='number', rho=1.):
        """Integrate the number of particles, surface area, volume, or mass
        between two diameters.

        Diameters should be in microns.
        """
        # First, let's find all bins within the range including the fringe bins!
        bin_idx = np.where((self.bins[:, -1] > dmin) & (self.bins[:,0] < dmax))

        _bins = self.bins[bin_idx]

        # generate an array of 'factors' to multiply by
        _f = np.ones(bin_idx[0].shape[0])

        for _i, _b in enumerate(_bins):
            f = 1.
            if _b[0] >= dmin and _b[-1] > dmax:
                _f[_i] = (dmax - _b[0]) / (_b[-1] - _b[0])
            elif _b[0] < dmin and _b[-1] >= dmin:
                _f[_i] = 1 - abs(_b[0] - dmin) / (_b[-1] - _b[0])

        # Get the data for the correct weight and multiply each row by the factors we just computed
        if weight == 'number':
            res = self.dn
        elif weight == 'surface':
            res = self.ds
        elif weight == 'volume':
            res = self.dv
        elif weight == 'mass':
            res = self.dv * rho

        # Grab the right data
        res = res.ix[:, bin_idx[0]].mul(_f).sum(axis=1)

        return res

    def heatmap(self):
        return heatmap(self.raw.index.values, self.midpoints, self.dndlogdp.T.values)

    def slice(self, start=None, end=None, inplace=False):
        """Slice the data between the start and end dates
        """
        if inplace:
            self.raw = self.raw[start:end]

            return True
        else:
            _tmp = _tmp.copy()

            _tmp.raw = _tmp.raw[start:end]

        return _tmp

    def resample(self, rs, inplace=False):
        """Resample the raw data"""
        obj_cols = self.raw.select_dtypes(include=['object']).resample(rs).first()
        num_cols = self.raw.select_dtypes(exclude=['object']).resample(rs).mean()

        # re-merge the two dataframes
        merged = pd.merge(num_cols, obj_cols, left_index=True, right_index=True, how='outer')

        if inplace:
            self.raw = merged

            return True
        else:
            _tmp = self.copy()
            _tmp.raw = merged

        return _tmp


def load_file(fpath, column=True, delimiter=',', encoding='ISO-8859-1', **kwargs):
    """Load an SMPS.dat file as exported using the TSI GUI"""
    # Get the number of lines to parse as meta info
    assert (delimiter in [',', '\t']), "Delimiter must be either tab or comma."

    _metacount = _get_linecount(fpath=fpath, keyword='Sample #', delimiter=delimiter)

    # Read in the meta data as a dictionary
    meta = pd.read_table(
                fpath,
                nrows=_metacount,
                delimiter=delimiter,
                header=None,
                encoding=encoding,
                error_bad_lines=False,
                warn_bad_lines=False,
                index_col=0).T.iloc[0,:].to_dict()

    unit = meta['Units']
    weight = meta['Weight']
    multiplier = float(meta['Channels/Decade'])

    if column is True:
        ts = pd.read_table(
                fpath,
                skiprows=_metacount,
                nrows=3,
                delimiter=delimiter,
                header=None,
                encoding=encoding).iloc[1:, 1:].T

        ts.columns = ['Date', 'Time']
        ts['timestamp'] = ts.apply(lambda x: pd.to_datetime("{} {}".format(x['Date'], x['Time'])), axis=1)

        # Retrieve the number of bins in the file
        nbins = _get_bin_count(fpath=fpath, delimiter=delimiter, encoding=encoding)

        # Read the table of raw data
        data = pd.read_table(
                    fpath,
                    skiprows=_metacount+4,
                    nrows=nbins,
                    delimiter=delimiter,
                    header=None,
                    encoding=encoding).astype(float)

        midpoints = data.iloc[:, 0].values
        data = data.iloc[:, 1:].T

        # Alter the column names
        data.columns = ['bin{}'.format(i) for i in range(nbins)]
        data.index = ts['timestamp']

        # Read the table of stats
        stats = pd.read_table(
                    fpath,
                    skiprows=(_metacount+4+nbins),
                    delimiter=delimiter,
                    header=None,
                    error_bad_lines=False,
                    warn_bad_lines=True,
                    encoding=encoding).T

        # Rename the columns to the first row
        stats = stats.rename(columns=stats.iloc[0])
        stats = stats.rename(columns=RENAMED_COLUMNS)

        # Drop the first row
        stats = stats.iloc[1:,:]

        # Set the index to timestamp
        stats.index = ts['timestamp']

        # Force the stats columns to be floats where needed
        for column in stats.columns:
            try:
                stats[column] = stats[column].astype(float)
            except: pass

        hist_cols = data.columns
        df = pd.merge(data, stats, left_index=True, right_index=True, how='outer')
    else:
        df = pd.read_table(
                fpath,
                skiprows=_metacount,
                delimiter=delimiter,
                encoding=encoding).dropna(how='all', axis=1)

        df.index = df.apply(lambda x: pd.to_datetime("{} {}".format(x['Date'], x['Start Time'])), axis=1)

        midpoints = []
        for col in df.columns:
            try:
                float(col)
                midpoints.append(col)
            except ValueError:
                pass

        hist_cols = ["bin{}".format(i) for i in range(len(midpoints))]

        df.rename(columns=dict(zip(midpoints, hist_cols)), inplace=True)

        midpoints = np.array([float(i) for i in midpoints])

        # Rename columns to be uniform and easier to parse
        df.rename(columns=RENAMED_COLUMNS, inplace=True)

    # Calculate bins
    bins = np.empty([len(midpoints), 3])
    bins.fill(np.NaN)

    # Make this shit more robust!
    bins[0, 0] = float(df['Lower Size'][0])
    bins[-1, -1] = float(df['Upper Size'][0])
    bins[:, 1] = midpoints

    for i in range(bins.shape[0] - 1):
        bins[i, 2] = round(math.pow(10, np.log10(bins[i, 0]) + (1./multiplier)), 4)
        bins[i + 1, 0] = bins[i, 2]

    return SMPS(data=df, meta=meta, bins=bins, bin_labels=hist_cols)

def load_sample(label="boston"):
    """Load a sample data file directly from GitHub.

    Options include: ['boston', 'chamber']

    """
    files = {
        'boston': {
            'uri': "https://raw.githubusercontent.com/dhhagan/py-smps/master/sample-data/boston_wintertime.txt",
            'column': False
        },
        'chamber': {
            'uri': "https://raw.githubusercontent.com/dhhagan/py-smps/master/sample-data/mit_chamber_sample_column.txt",
            'column': True
        }
    }

    return load_file(files[label]['uri'], column=files[label]['column'])
