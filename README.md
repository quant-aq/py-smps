[![Code Climate](https://codeclimate.com/github/dhhagan/py-smps/badges/gpa.svg)](https://codeclimate.com/github/dhhagan/py-smps)
[![Test Coverage](https://codeclimate.com/github/dhhagan/py-smps/badges/coverage.svg)](https://codeclimate.com/github/dhhagan/py-smps/coverage)
[![Issue Count](https://codeclimate.com/github/dhhagan/py-smps/badges/issue_count.svg)](https://codeclimate.com/github/dhhagan/py-smps)

# py-smps
Python library for the analysis and visualization of data from a Scanning Mobility Particle Sizer (SMPS) and other similar instruments (SEMS, OPC's).

## Dependencies

  * pandas
  * numpy
  * seaborn

## Installation

To install directly from GitHub using pip:

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

I will eventually get to this, but for now you can check out the examples [here](/examples).

## Contributing to Development

I threw this together because I was tired of analyzing SMPS data with proprietary software. Please help me contribute?

If there is a feature you would like to see or a bug you would like to report, please open an issue. I will try to get to things as promptly as possible. Otherwise, feel free to send PR's!


## Colorbar Information

  * [matplotlib colorbars](http://matplotlib.org/examples/color/colormaps_reference.html)
  * [seaborn color palette](http://seaborn.pydata.org/tutorial/color_palettes.html)
