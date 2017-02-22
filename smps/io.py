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
    """"""
    def __init__(self, df, meta, bins, bin_labels):
        self.df = df
        self.meta = meta
        self.bins = bins
        self.midpoints = bins[:, 1]
        self.bin_labels = bin_labels

        self.histogram = self.df[self.bin_labels]

    def histogram(self):
        return self.df[self.bin_labels]

    def heatmap(self):
        return heatmap(self.df.index.values, self.midpoints, self.histogram.T.values)

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

    return SMPS(df, meta, bins, hist_cols)
