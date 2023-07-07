# Getting Started

New to _py-smps_? Don't worry, you've found the perfect place to get started.

## Installation

_py-smps_ can be installed by running:

```sh
$ pip install py-smps
```

It requires Python 3.7+ to run. 

If you can't wait for the latest _hotness_ and want to install from GitHub, use:

```sh
$ pip install git+https://github.com/quant-aq/py-smps
```

## Basic Usage

To get started, it's as simple as loading the library and reading in some data (or using one of the example data sets).

```python
import smps

obj = smps.io.load_sample("boston")

# Compute the stats
stats = obj.stats(weight='number')
```



## Next Steps

Take a look at the examples and guides in the **Guides** section for an in-depth overview of what the _py-smps_ library is capable of. 