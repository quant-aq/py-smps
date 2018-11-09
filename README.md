[![PyPI version](https://badge.fury.io/py/smps.svg)](https://badge.fury.io/py/smps)
[![Build Status](https://travis-ci.org/dhhagan/py-smps.svg?branch=master)](https://travis-ci.org/dhhagan/py-smps)
[![Coverage Status](https://coveralls.io/repos/github/dhhagan/py-smps/badge.svg?branch=master)](https://coveralls.io/github/dhhagan/py-smps?branch=master)


# py-smps
Python library for the analysis and visualization of data from a Scanning Mobility Particle Sizer (SMPS) and other similar instruments (SEMS, OPC's).

If you use this for scientific research, please cite? It's all I've got...

## Dependencies

  * pandas
  * numpy
  * scipy
  * seaborn
  * statsmodels

## Python Versions

Currently, only Python3.3+ is supported. If you are using/relying on Python2.7, it is time to switch!

## Installation

To install from PyPi:

    $ pip install py-smps [--upgrade]

To install the edge release directly from GitHub:

    pip install git+https://github.com/dhhagan/py-smps.git

Feel free to download or clone the repository as well and install from source.

    python3 setup.py install (--upgrade)

## Unittests

Unittests can be run by issuing the following command from within the main repo:

    $ python3 setup.py test

They can also be run with coverage (if installed) by running the following in succession:

    $ coverage run --source smps setup.py test
    $ coverage report -m


## Documentation

I will eventually finish writing nice/new documentation, but for now you can check out the examples [here](/examples) which cover most use cases.

## Contributing to Development

I threw this together because I was tired of analyzing SMPS data with proprietary software. Please help me contribute?

If there is a feature you would like to see or a bug you would like to report, please open an issue. I will try to get to things as promptly as possible. Otherwise, feel free to send PR's!


## Colorbar Information

  * [matplotlib colorbars](http://matplotlib.org/examples/color/colormaps_reference.html)
  * [seaborn color palette](http://seaborn.pydata.org/tutorial/color_palettes.html)
