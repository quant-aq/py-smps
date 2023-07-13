#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import pandas as pd
import numpy as np
import math
import copy
import joblib

from .utils import _get_bin_count, _get_linecount, make_bins
from .plots import heatmap
from .models import SMPS


__all__ = ["smps_from_txt", "load_sample"]


def smps_from_txt(fpath, column=True, delimiter=',', as_dict=False, **kwargs):
    """Create an SMPS object directory from a text file via TSIs AIM software.
    
    Parameters
    ----------
    fpath : str
            The filepath where the file to read is located.
    column : bool, default=True
            If your data was exported in 'column' format, leave as True. If exported
            in 'row' format, set to False.
    delimiter : str, default=','
            Sets the delimiter to expect when reading the data file.
    as_dict : bool, default=False
            Set to true to return a dict, otherwise returns an SMPS object.
    
    Returns
    -------
    smps.models.SMPS or dict
    
    See Also
    --------
    smps.models.SMPS
    
    Examples
    --------
    
    Load a data file with all defaults:
    
    >>> obj = smps.io.smps_from_text(filepath)
    
    Load a similar file, but exported in row format:
    
    >>> obj = smps.io.smps_from_text(filepath, column=False)
    
    """
    assert(delimiter in [',', '\t']), "The delimiter must \
    be either a comma or tab"

    encoding = kwargs.pop("encoding", "ISO-8859-1")

    # determine the number of rows of meta information; can be defined
    meta_num_lines = kwargs.pop("meta_num_lines", None)
    if not meta_num_lines:
        meta_num_lines = _get_linecount(
            fpath, 
            keyword='Sample #', 
            encoding=encoding, 
            delimiter=delimiter
        )

    # read and store the meta information
    meta = pd.read_table(
        fpath, 
        nrows=meta_num_lines, 
        delimiter=delimiter, 
        header=None, 
        encoding=encoding,
        on_bad_lines='warn',
        index_col=0
    ).T.iloc[0,:].to_dict()

    # grab the number of channels per decade - sets the 
    # multiplier for determining bin spacings later on
    mult = float(meta['Channels/Decade'])

    # determine the weight and units
    units = meta['Units'].lower()
    weight = meta['Weight'].lower()

    # read the data
    if column:
        ts = pd.read_table(
            fpath,
            skiprows=meta_num_lines,
            nrows=3,
            delimiter=delimiter,
            header=None,
            encoding=encoding).iloc[1:, 1:].T

        ts.columns = ['Date', 'Time']

        # set the timestamp
        ts['timestamp'] = ts.apply(
            lambda x: pd.to_datetime(
                "{} {}".format(x['Date'], 
                               x['Time']
                              )
            ), axis=1
        )

        # Retrieve the number of bins in the file
        nbins = _get_bin_count(
            fpath=fpath, 
            delimiter=delimiter, 
            encoding=encoding
        )

        # Read the table of raw data
        data = pd.read_table(
                    fpath,
                    skiprows=meta_num_lines+4,
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
                    skiprows=(meta_num_lines+4+nbins),
                    delimiter=delimiter,
                    header=None,
                    on_bad_lines='warn',
                    encoding=encoding).T

        # Rename the columns to the first row
        stats = stats.rename(columns=stats.iloc[0])

        # Drop the first row
        stats = stats.iloc[1:,:]

        # Set the index to timestamp
        stats.index = ts['timestamp']

        # Force the stats columns to be floats where needed
        for column in stats.columns:
            try:
                stats[column] = stats[column].astype(float)
            except: pass

        bin_labels = list(data.columns)

        data = pd.merge(
            data, 
            stats, 
            left_index=True, 
            right_index=True, 
            how='outer'
        )
    else:
        data = pd.read_table(
                    fpath,
                    skiprows=meta_num_lines,
                    delimiter=delimiter,
                    encoding=encoding).dropna(how='all', axis=1)

        data.index = data.apply(
            lambda x: pd.to_datetime("{} {}".format(
                x['Date'], 
                x['Start Time']
            )), axis=1
        )

        # delete some un-needed columns
        del data['Date'], data['Start Time']

        # grab the midpoint diameters
        midpoints = []
        for col in data.columns:
            try:
                float(col)
                midpoints.append(col)
            except ValueError:
                pass

        bin_labels = ["bin{}".format(i) for i in range(len(midpoints))]

        # rename the bins
        data.rename(
            columns=dict(zip(midpoints, bin_labels)), 
            inplace=True
        )

        midpoints = np.array([float(i) for i in midpoints])

    # generate the bins
    low_col_name = kwargs.pop(
        "lower_column_label", 
        [c for c in data.columns if "lower" in c.lower()][0]
    )
    upper_col_name = kwargs.pop(
        "upper_column_label", 
        [c for c in data.columns if "upper" in c.lower()][0]
    )
    bound_left = float(data[low_col_name][0])
    bound_right = float(data[upper_col_name][0])

    meta['Lower Size (nm)'], meta['Upper Size (nm)'] \
    = bound_left, bound_right

    # make the bin array
    bins = make_bins(
        midpoints=midpoints, 
        lb=bound_left, 
        ub=bound_right,
        channels_per_decade=int(meta['Channels/Decade'])
    )

    # rename the index
    data.index.rename("timestamp", inplace=True)

    if as_dict:
        return dict(
                meta=meta,
                units=units,
                weight=weight,
                data=data,
                bins=bins,
                bin_labels=bin_labels,
                bin_prefix='bin')
    else:
        return SMPS(data=data, meta=meta, bins=bins, 
                    bin_labels=bin_labels, units=units, weight=weight)


def load_sample(label="boston"):
    """Load one of the example datasets provided.
    
    Parameters
    ----------
    label : str, default='boston'
            Choose the example dataset to load. Should be one of ['boston', 'chamber']
    
    Returns
    -------
    smps.models.SMPS
    
    See Also
    --------
    smps.models.SMPS
    
    Examples
    --------
    
    Load the boston example dataset:
    
    >>> obj = smps.io.load_sample('boston')
    
    """
    assert(label in ["boston", "chamber"]), "Invalid option chosen \
    for the label"

    files = {
        'boston': {
            'uri': ("https://raw.githubusercontent.com/dhhagan/py-"
            + "smps/master/sample-data/boston_wintertime.txt"),
            'column': False
        },
        'chamber': {
            'uri': ("https://raw.githubusercontent.com/dhhagan/py-"
            + "smps/master/sample-data/mit_chamber_sample_column.txt"),
            'column': True
        }
    }

    m = smps_from_txt(fpath=files[label]['uri'], column=files[label]['column'], as_dict=True)

    # convert to an SMPS instance
    return SMPS(
        data=m['data'], 
        bins=m['bins'], 
        meta=m['meta'], 
        bin_labels=m['bin_labels'], 
        weight=m['weight'], 
        units=m['units']
    )
