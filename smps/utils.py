#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import requests

def _get_bin_count(fpath, encoding='ISO-8859-1'):
    """Get the number of bins in the file."""
    bins = 0
    if 'http' in fpath:
        req = requests.get(fpath)

        for line in req.iter_lines():
            try:
                if float(line.decode(encoding).split(',')[0]):
                    bins += 1
            except: pass
    else:
        with open(fpath, 'r', encoding=encoding) as f:
            for line in f:
                try:
                    if float(line.split(',')[0]):
                        bins += 1
                except: pass

    return bins

def _get_linecount(fpath, keyword, delimiter=',', encoding='ISO-8859-1'):
    """Return the line number in a file where the first item is `keyword`.
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
    return x if x % 100 == 0 else x + 100 - x % 100

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
