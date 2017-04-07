"""
"""
import pandas as pd
import numpy as np
import math

from .utils import _get_bin_count
from .plots import heatmap

SMPS_STATS_COLUMN_NAMES = [
    'Scan Up Time',
    'Retrace Time',
    'Down Scan First',
    'Scans Per Sample',
    'Impactor Type',
    'Sheath Flow',
    'Aerosol Flow',
    'CPC Inlet Flow',
    'CPC Sample Flow',
    'Low Voltage',
    'High Voltage',
    'Lower Size',
    'Upper Size',
    'Density',
    'Title',
    'Status Flag',
    'td',
    'tf',
    'D50',
    'Median',
    'Mean',
    'Geo. Mean',
    'Mode',
    'Geo Std Dev',
    'Total Concentration',
    'Comment'
]

class SMPS(object):
    """Assumes data is always fed as dNdlogDp"""
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

    @property
    def dlogdp(self):
        """Return the bin midpoints"""
        return np.log10(self.bins[:, -1]) - np.log10(self.bins[:, 0])

    @property
    def dndlogdp(self):
        """"""
        return self.raw[self.bin_labels]

    @property
    def dn(self):
        return self.dndlogdp * self.dlogdp

    @property
    def ds(self):
        return self.dsdlogdp.mul(self.dlogdp)

    @property
    def dv(self):
        """Returns the volume in each bin"""
        return self.dvdlogdp.mul(self.dlogdp)

    @property
    def dsdlogdp(self):
        """Return dSdlogDp"""
        return self.dndlogdp.mul(self.s_multiplier)

    @property
    def dvdlogdp(self):
        """Return dVdlogDp"""
        return self.dndlogdp.mul(self.v_multiplier)

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

def load_file(fpath, column=True, **kwargs):
    """"""
    delim    = kwargs.pop('delimiter', ',')
    encoding = kwargs.pop('encoding', 'ISO-8859-1')

    # Read in the meta data as a dictionary
    meta = pd.read_table(
                fpath,
                nrows=15,
                delimiter=delim,
                header=None,
                encoding=encoding,
                index_col=0).T.ix[1,:].to_dict()

    unit = meta['Units']
    weight = meta['Weight']
    multiplier = float(meta['Channels/Decade'])

    if column is True:
        ts = pd.read_table(
                fpath,
                skiprows=15,
                nrows=3,
                delimiter=delim,
                header=None,
                encoding=encoding).ix[:, 1:].T

        ts.columns = ['Sample #', 'Date', 'Start Time']
        ts['timestamp'] = ts.apply(lambda x: pd.to_datetime("{} {}".format(x['Date'], x['Start Time'])), axis=1)

        nbins = _get_bin_count(fpath, encoding)

        data = pd.read_table(
                    fpath,
                    skiprows=19,
                    nrows=nbins,
                    delimiter=delim,
                    header=None,
                    encoding=encoding).astype(float)

        midpoints = data.ix[:, 0].values
        data = data.ix[:, 1:].T

        # Alter the column names
        data.columns = ['bin{}'.format(i) for i in range(nbins)]
        data.index = ts['timestamp']

        stats = pd.read_table(
                    fpath,
                    skiprows=(19+nbins),
                    delimiter=delim,
                    header=None,
                    encoding=encoding).ix[:, 1:].T

        stats.columns = SMPS_STATS_COLUMN_NAMES
        stats.index = ts['timestamp']

        hist_cols = data.columns
        df = pd.merge(data, stats, left_index=True, right_index=True, how='outer')
    else:
        df = pd.read_table(
                fpath,
                skiprows=15,
                delimiter=delim,
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

    # Rename some columns
    df.rename(columns={'Total Conc.(#/cm³)': 'Total Conc.'}, inplace=True)

    # Calculate bins
    bins = np.empty([len(midpoints), 3])
    bins.fill(np.NaN)

    bins[0, 0] = float(df['Lower Size(nm)'][0])
    bins[-1, -1] = float(df['Upper Size(nm)'][0])
    bins[:, 1] = midpoints

    for i in range(bins.shape[0] - 1):
        bins[i, 2] = round(math.pow(10, np.log10(bins[i, 0]) + (1./multiplier)), 4)
        bins[i + 1, 0] = bins[i, 2]

    return SMPS(data=df, meta=meta, bins=bins, bin_labels=hist_cols)

def load_sample(label="boston"):
    """Load sample data.

    ** THIS DOESN'T WORK **
    """
    f = load_file("../sample-data/boston_wintertime.txt", column=False)

    return f
