"""
"""
import pandas as pd

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
    def __init__(self, df, meta, midpoints, bin_labels):
        self.df = df
        self.meta = meta
        self.midpoints = midpoints
        self.bin_labels = bin_labels

        self.histogram = self.df[self.bin_labels]

    def histogram(self):
        return self.df[self.bin_labels]

    def heatmap(self):
        return heatmap(self.df.index.values, self.midpoints, self.histogram.T.values)

def load_smps_file(fpath, **kwargs):
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

    # Read in the data
    data = pd.read_table(
            fpath,
            #nrows=3,
            skiprows=15,
            delimiter=delim,
            header=None,
            encoding=encoding)

    ts = data.ix[1:2, 1:].T
    ts.columns = ['Date', 'Time']

    # Turn into a timeseries that will become the index
    ts = ts.apply(lambda x: pd.to_datetime("{} {}".format(x['Date'], x['Time'])), axis=1)

    # Get the histogram
    nbins = _get_bin_count(fpath, encoding)
    midpoints = data.ix[4:4 + nbins - 1, 0].astype(float).values

    hist = data.ix[4:4 + nbins - 1, 1:].T
    hist.index = ts
    hist.columns = ["bin{}".format(i) for i in range(nbins)]
    hist = hist.astype(float)

    stats = data.ix[4 + nbins:, 1:].T
    stats.index = ts
    stats.columns = SMPS_STATS_COLUMN_NAMES

    df = pd.merge(hist, stats, left_index=True, right_index=True, how='outer')

    return SMPS(df, meta, midpoints, hist.columns)
