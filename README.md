[![PyPI version](https://badge.fury.io/py/py-smps.svg)](https://badge.fury.io/py/py-smps)
[![Build Status](https://travis-ci.org/dhhagan/py-smps.svg?branch=master)](https://travis-ci.org/dhhagan/py-smps)
[![Coverage Status](https://coveralls.io/repos/github/dhhagan/py-smps/badge.svg?branch=master)](https://coveralls.io/github/dhhagan/py-smps?branch=master)


# py-smps
Python library for the analysis and visualization of data from a Scanning Mobility Particle Sizer (SMPS) and other particle sizing instruments (SEMS, OPC's).

## Dependencies

The full list of dependencies can be found in the `pyproject.toml` file and are summarized below:

```py
python = ">=3.8, <3.12"
statsmodels = "^0.13"
seaborn = "^0.10"
joblib = "^1.0"
requests = "^2.24"
scipy = "^1.9"
numpy = "^1.23.2"
pandas = "^1.4"
```

As of `v1.2.0a0`, the library should be compatible with Apple silicone (tested on both M1 and M2).

## Python Versions

Python3.8 through Python3.11 are currently supported.

## Installation

To install from PyPi:

    $ pip install py-smps [--upgrade]

If you'd like the latest pre-release:

    $ pip install py-smps --pre [--upgrade]

To install the edge release directly from GitHub:

    pip install git+https://github.com/quant-aq/py-smps.git

## Unittests

Unittests can be run by issuing the following command from within the main repo:

```sh
$ poetry run pytest -s tests/ --ignore=tests/datafiles
```


## Documentation

Documentation is available [here](https://quant-aq.github.io/py-smps/). Docs are built using Sphinx and can be built locally by doing the following:

```sh
$ cd docs/
$ make clean
$ make guides
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

