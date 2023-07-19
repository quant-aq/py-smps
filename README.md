[![PyPI version](https://badge.fury.io/py/py-smps.svg)](https://badge.fury.io/py/py-smps)
<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-1-orange.svg?style=flat-square)](#contributors-)
<!-- ALL-CONTRIBUTORS-BADGE:END -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Coverage Status](https://coveralls.io/repos/github/dhhagan/py-smps/badge.svg?branch=master)](https://coveralls.io/github/dhhagan/py-smps?branch=master)
[![ci.tests](https://github.com/quant-aq/py-smps/actions/workflows/test-and-report.yml/badge.svg)](https://github.com/quant-aq/py-smps/actions/workflows/test-and-report.yml)


# py-smps

py-smps is a Python data analysis library built for analyzing size-resolved aerosol data from a variety of aerosol sizing instruments (e.g., Scanning Mobility Particle Sizer, Optical Particle Counters).


**NOTE: As of `v1.2.0`, the library is compatible with Apple silicone (M1, M2 chips).**

# Installation

Official releases of `py-smps` can be installed from [PyPI](https://pypi.org/project/py-smps/):

    $ pip install py-smps [--upgrade]

If you'd like the latest pre-release:

    $ pip install py-smps --pre [--upgrade]

To install the edge release directly from GitHub:

    pip install git+https://github.com/quant-aq/py-smps.git

# Dependencies

## Supported Python versions
- Python 3.8+

## Mandatory Dependencies

The full list of dependencies can be found in the [`pyproject.toml`](pyproject.toml) file.

# Development

## Testing

Tests can be run by issuing the following command from within the main repo:

```sh
$ poetry run pytest -s tests/ --ignore=tests/datafiles
```

## Contributing to Development

We welcome all contributions from the community in the form of issues reporting, feature requests, bug fixes, etc.

If there is a feature you would like to see or a bug you would like to report, please open an issue. We will try to get to things as promptly as possible. Otherwise, feel free to send PR's!

### Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://www.davidhagan.me"><img src="https://avatars.githubusercontent.com/u/4370504?v=4?s=100" width="100px;" alt="David H Hagan"/><br /><sub><b>David H Hagan</b></sub></a><br /><a href="https://github.com/quant-aq/py-smps/commits?author=dhhagan" title="Code">💻</a> <a href="https://github.com/quant-aq/py-smps/commits?author=dhhagan" title="Tests">⚠️</a> <a href="https://github.com/quant-aq/py-smps/commits?author=dhhagan" title="Documentation">📖</a> <a href="#content-dhhagan" title="Content">🖋</a></td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->


# Documentation

Documentation is available [here](https://quant-aq.github.io/py-smps/). To build locally, you must first install [pandoc](https://pandoc.org/). Docs are built using Sphinx and can be built locally by doing the following:

```sh
# Activate the virtualenv
$ poetry shell

# Build the docs
$ cd docs/
$ make clean
$ make html
$ cd ..
```

Then, you can navigate to your local directory at `docs/build/html/` and open up the `index.html` file in your preferred browser window.

## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!