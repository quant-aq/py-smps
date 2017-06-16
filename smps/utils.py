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
                if float(line.decode('utf-8').split(',')[0]):
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

def roundup(x):
    return x if x % 100 == 0 else x + 100 - x % 100
