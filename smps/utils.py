#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import requests

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


RENAMED_COLUMNS = {
    'Scan Up Time(s)': 'Scan Up Time',
    'Retrace Time(s)': 'Retrace Time',
    'Impactor Type(cm)': 'Impactor Type',
    'Sheath Flow(lpm)': 'Sheath Flow',
    'Aerosol Flow(lpm)': 'Aerosol Flow',
    'CPC Inlet FLow(lpm)': 'CPC Inlet Flow',
    'CPC Sample Flow(lpm)':'CPC Sample Flow',
    'Lower Size(nm)': 'Lower Size',
    'Upper Size(nm)': 'Upper Size',
    'Density(g/cc)': 'Density',
    'td(s)': 'td',
    'tf(s)': 'tf',
    'D50(nm)': 'D50',
    'Median(nm)': 'Median',
    'Mean(nm)': 'Mean',
    'Median(nm)': 'Median',
    'Geo. Mean(nm)': 'GM',
    'Mode(nm)': 'Mode',
    'Geo. Std. Dev.': 'GSD',
    'Total Conc.(#/cm³)': 'Total Conc.',
    'Total Concentration': 'Total Conc.'
 }

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
    'GM',
    'Mode',
    'GSD',
    'Total Conc.',
    'Comment'
    ]
