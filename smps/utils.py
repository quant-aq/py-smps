#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import requests
import math

def make_bins(nbins, bound_left=None, bound_right=None, midpoints=None,
                channels_per_decade=None, mean='gm'):
    """

    :param mean: ['gm', 'am']

    Example:

    >>> make_bins(
    >>>  len(obj['midpoints']),
    >>>  midpoints=obj['midpoints'],
    >>>  bound_left=obj['bound_left'],
    >>>  bound_right=obj['bound_right'],
    >>>  channels_per_decade=64)
    """
    # initialize the bins array which will contain an nbinsx3 matrix
    bins = np.empty((nbins, 3))
    bins.fill(np.NaN)

    # if midpoints is an array, set it
    if type(midpoints) in (np.ndarray, list):
        bins[:, 1] = midpoints

    # set the left bounds
    if type(bound_left) in (np.ndarray, list):
        bins[:, 0] = bound_left
    else:
        bins[0, 0] = bound_left

    # set the right bounds
    if type(bound_right) in (np.ndarray, list):
        bins[:, -1] = bound_right
    else:
        bins[-1, -1] = bound_right

    # iterate over all bins and set the missing data
    for i in range(bins.shape[0] - 1):
        if bins[i, 1] == np.nan: # set the midpoint as either the GM or the AM of the two
            if mean == 'am':
                bins[i, 1] = np.mean([bins[i, 0], bins[i, -1]])
            else:
                bins[i, 1] = 1
        else:
            bins[i, -1] = round(math.pow(10, np.log10(bins[i, 0]) + 1./channels_per_decade), 4)
            bins[i+1, 0] = bins[i, -1]

    return bins

def _get_bin_count(fpath, delimiter=',', encoding='ISO-8859-1'):
    """Get the number of bins in the file."""
    bins = 0
    if 'http' in fpath:
        req = requests.get(fpath)

        for line in req.iter_lines():
            try:
                if float(line.decode(encoding).split(delimiter)[0]):
                    bins += 1
            except: pass
    else:
        with open(fpath, 'r', encoding=encoding) as f:
            for line in f:
                try:
                    if float(line.split(delimiter)[0]):
                        bins += 1
                except: pass

    return bins

def _get_linecount(fpath, keyword, delimiter=',', encoding='ISO-8859-1'):
    """Return the line number in a file where the first item is `keyword`"""
    linecount = 0

    if 'http' in fpath:
        req = requests.get(fpath)

        for line in req.iter_lines():
            startswith = line.decode(encoding).split(delimiter)[0]

            if startswith == keyword:
                break

            linecount += 1
    else:
        with open(fpath, 'r', encoding=encoding) as f:
            for line in f:
                startswith = line.split(delimiter)[0]

                if startswith == keyword:
                    break

                linecount += 1

    return linecount

def roundup(x):
    return x if x % 100 == 0 else x + 100 - x % 100




class Table(object):
    def __init__(self, max_width=80, **kwargs):
        self.text = ""
        self.max_width_chars = max_width

    def add_title(self, title):
        self.text = title.center(self.max_width_chars) + self.text

        return

    def add_border(self, char='='):
        self.text += '\n' + char*self.max_width_chars

        return

    def _center_text(self, text, width=None):
        if width is None:
            width = int(self.max_width_chars / 4)

        return text.center(width)

    def add_header(self):
        self.text += "\n" + self._center_text("") + self._center_text("N (cm-3)")
        self.text += self._center_text("GM (nm)") + self._center_text("GSD")

        return

    def add_row(self, label, fields, errors):
        self.text += "\n" + self._center_text(label)

        field1 = "{:.2e} ({:.1e})".format(fields[0], errors[0])
        field2 = "{:.2f} ({:.1e})".format(fields[1], errors[1])
        field3 = "{:.2f} ({:.1e})".format(fields[2], errors[2])

        self.text += self._center_text(field1) + self._center_text(field2) + \
            self._center_text(field3)

        return

    def __repr__(self):
        return self.text
