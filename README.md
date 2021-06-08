[![PyPI version](https://badge.fury.io/py/py-smps.svg)](https://badge.fury.io/py/py-smps)
[![Build Status](https://travis-ci.org/dhhagan/py-smps.svg?branch=master)](https://travis-ci.org/dhhagan/py-smps)
[![Coverage Status](https://coveralls.io/repos/github/dhhagan/py-smps/badge.svg?branch=master)](https://coveralls.io/github/dhhagan/py-smps?branch=master)


# py-smps
Python library for the analysis and visualization of data from a Scanning Mobility Particle Sizer (SMPS) and other particle sizing instruments (SEMS, OPC's).

## Dependencies

  * pandas
  * numpy
  * scipy
  * seaborn
  * statsmodels

## Python Versions

Currently, Python3.7+ is supported and we test against Python 3.7, 3.8, and 3.9.

## Installation

To install from PyPi:

    $ pip install py-smps [--upgrade]

To install the edge release directly from GitHub:

    pip install git+https://github.com/quant-aq/py-smps.git

## Unittests

Unittests can be run by issuing the following command from within the main repo:

```sh
$ poetry run pytest tests/ --ignore=tests/datafiles
```


## Documentation

Documentation is available [here](https://quant-aq.github.io/py-smps/). Docs are built using Sphinx and can be built locally by doing the following:

```sh
$ cd docs/
$ make clean
$ make html
$ cd ..
```

Then, you can navigate to your local directory at `docs/_build/html/` and open up the `index.html` file in your preferred browser window.


## Contributing to Development

We welcome all contributions from the community in the form of issues reporting, feature requests, bug fixes, etc.

If there is a feature you would like to see or a bug you would like to report, please open an issue. We will try to get to things as promptly as possible. Otherwise, feel free to send PR's!


## Colorbar Information

  * [matplotlib colorbars](http://matplotlib.org/examples/color/colormaps_reference.html)
  * [seaborn color palette](http://seaborn.pydata.org/tutorial/color_palettes.html)
