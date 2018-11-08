#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import requests
import math
from scipy.stats.mstats import gmean

def make_bins(**kwargs):
    """
    """
    boundaries = kwargs.pop("boundaries", None)
    boundaries_left = kwargs.pop("boundaries_left", None)
    boundaries_right = kwargs.pop("boundaries_right", None)
    lb = kwargs.pop("lb", None)
    ub = kwargs.pop("ub", None)
    midpoints = kwargs.pop("midpoints", None)
    cpd = kwargs.pop("channels_per_decade", 64)
    mean_calc = kwargs.pop("mean_calc", "am")

    # initialize bins
    bins = None

    if midpoints is not None:
        if lb is None or ub is None:
            raise Exception("A lower and upper bound must be set")

        bins = np.empty((midpoints.shape[0], 3))
        bins.fill(np.NaN)

        # fill the midpoints
        bins[:, 1] = midpoints

        # fill the bounds
        bins[0, 0] = lb
        bins[-1, -1] = ub

        # iterate and calculate the bounds
        for i in range(bins.shape[0] - 1):
            bins[i, 2] = round(math.pow(10, np.log10(bins[i, 0]) + 1./cpd), 4)
            bins[i+1, 0] = bins[i, 2]

        return bins

    if boundaries is not None:
        bins = np.empty((boundaries.shape[0]-1, 3))
        bins.fill(np.NaN)

        bins[:, 0] = boundaries[0:-1]
        bins[:, 2] = boundaries[1:]

    elif boundaries_left is not None:
        if boundaries_right is None:
            raise Exception("Missing attribute: `boundaries_right`")

        assert(boundaries_left.shape[0] == boundaries_right.shape[0]), \
            "Boundaries must be the same dimensions."

        bins = np.empty((boundaries_left.shape[0], 3))
        bins.fill(np.NaN)

        bins[:, 0] = boundaries_left
        bins[:, 2] = boundaries_right

    else:
        raise Exception("Not enough information to compute.")

    # calculate the midpoints
    assert(mean_calc in ("gm", "am")), "Invalid mean calculation method."

    if mean_calc == 'am':
        bins[:, 1] = (bins[:, 0] + bins[:, 2]) / 2
    else:
        bins[:, 1] = [gmean([x, y]) for x, y in zip(bins[:, 0], bins[:, 2])]

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
