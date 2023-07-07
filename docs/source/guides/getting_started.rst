Getting Started
===============

Are you new to `py-smps`? Don't worry, you've found the perfect place to get started!

Installation
------------

The official release of `py-smps` can be installed from PyPI:

.. code-block:: console

    $ pip install py-smps [--upgrade]

If you'd like to install the most recent pre-release:

.. code-block:: console

    $ pip install py-smps --pre [--upgrade]

To install the most up-to-date development version directly from GitHub:

.. code-block:: console

    $ pip install git+https://github.com/quant-aq/py-smps.git@<tag-or-version-or-branch>


Basic Usage
-----------

To get up and running, it is as simple as loading the library and reading in some data (or using 
one of the example datasets available).

.. code-block:: python

    import smps

    # Load the sample Boston dataset
    obj = smps.io.load_sample("boston")

    # Compute the aerosol statistics, weighted by number
    stats = obj.stats(weight='number')


Next Steps
----------

Take a look at the examples and guides in the :ref:`Guides` section for an in-depth overview of what 
the `py-smps` library is capable of!