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
    'make_bins` computes the upper and lower bounds and
    midpoints of all the bins given either
    * `boundaries`
    * `boundaries_left` and `boundaries_right`
    * or `lb`, `ub`, and `midpoints` (`channels_per_decade` optional).
    
    :param boundaries: A list containing, for each bin, a list of 
        the lower and upper bound for the bin.
    :type boundaries: list of lists, or numpy array
    :param boundaries_left: A list of lower boundaries for each bin.
    :type boundaries_left: list
    :param boundaries_right: A list of upper boundaries for each 
        bin.
    :type boundaries_right: list
    :param lb: The lower bound of the lowermost bin.
    :type lb: float
    :param ub: The upper bound of the uppermost bin.
    :type ub: float
    :param midpoints: A list of midpoints for each bin.
    :type midpoints: list
    :param channels_per_decade: A measure of bin width on the log 
        scale.
    :type channels_per_decade: float
    :param mean_calc: Determins wether the geometric or arethmetic 
        mean should be used to compute the midpoints.
    :type mean_calc: {"gm", "am"}
    :return: For each bin a row is given with the lower boundary, 
        the midpoint, and the upper bounary for that bin.
    :rtype: numpy array
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
            bins[i, 2] = round(
                math.pow(10, np.log10(bins[i, 0]) + 1./cpd), 
                4
            )
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
        bins[:, 1] = [gmean([x, y]) for x, y in zip(
            bins[:, 0], 
            bins[:, 2]
        )]

    return bins

def _get_bin_count(fpath, delimiter=',', encoding='ISO-8859-1'):
    """
    Gets the number of bins in the file.
    
    :param fpath: A path or url for the file (if a url it must 
        include `http`, and if a file path it must not contain 
        `http`).
    :type fpath: string
    :param delimiter: The delimiter between items in the file, 
        defaults to ','.
    :type delimiter: string
    :param encoding: The encoding for the file, defaults to 
        'ISO-8859-1'.
    :type encoding: string
    :return: The number of bins in the file.
    :rtype: int
    """
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
    """
    Return the line number in a file where the first item is 
    `keyword`. If there is no such line, it will return the total 
    number of lines in the file.
    
    :param fpath: A path or url for the file (if a url it must 
        include `http`, and if a file path it must not 
        contain `http`).
    :type fpath: string
    :param keyword: The string to look for.
    :type keyword: string
    :param delimiter: The delimiter between items in the file, 
        defaults to ','.
    :type delimiter: string
    :param encoding: The encoding for the file, defaults to 
        'ISO-8859-1'.
    :type encoding: string
    :return: The line number in the file where the first item is 
        `keyword`.
    :rtype: int
    """
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
    """
    round up to a multiple of 100
    
    :param x: The number to be rounded up.
    :type x: float
    :return: `x` rounded up to a multiple of 100.
    :rtype: float
    """
    return x if x % 100 == 0 else x + 100 - x % 100




class Table(object):
    """
    A table which can be constructed and then printed into 
    the terminal.
    
    :param max_width: The maximum width, in characters, of the 
        table to be printed in the terminal, defaults to 80.
    :type max_width: int
    """
    def __init__(self, max_width=80, **kwargs):
        self.text = ""
        self.max_width_chars = max_width

    def add_title(self, title):
        """
        Add a title to the top of the table.
        
        :param title: The title to be added.
        :type title: string
        """
        self.text = title.center(self.max_width_chars) + self.text

        return

    def add_border(self, char='='):
        """
        Add a border to the bottom of the table to visually 
        delimit between data.
        
        :param char: The character to be used for the border.
        :type char: string
        """
        assert(len(char)==1), "`char` must be a single character"
        self.text += '\n' + char*self.max_width_chars

        return

    def _center_text(self, text, width=None):
        """
        Pad text on each side with spaces so that it takes up 
        `width` characters.
        
        :param text: The text to center.
        :type text: string
        :param width: The width, in characters, in which to center 
            the text, defaults to a quarter of `max_width`.
        :type width: int
        :return: The padded text.
        :rtype: string        
        """
        if width is None:
            width = int(self.max_width_chars / 4)

        return text.center(width)

    def add_header(self):
        """
        Appends a header to the bottom of the table with four 
        sections: blank, "N (cm-3)", "GM (nm)", and "GSD".
        """
        self.text += (
            "\n" 
            + self._center_text("") 
            + self._center_text("N (cm-3)")
        )
        self.text += (
            self._center_text("GM (nm)") 
            + self._center_text("GSD")
        )

        return

    def add_row(self, label, fields, errors):
        self.text += "\n" + self._center_text(label)

        field1 = "{:.2e} ({:.1e})".format(fields[0], errors[0])
        field2 = "{:.2f} ({:.1e})".format(fields[1], errors[1])
        field3 = "{:.2f} ({:.1e})".format(fields[2], errors[2])

        self.text += (
            self._center_text(field1) 
            + self._center_text(field2) 
            + self._center_text(field3)
        )

        return

    def __repr__(self):
        return self.text
