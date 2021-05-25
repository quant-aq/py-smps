Contribution
============

We welcome all contributions to `py-smps`.

Support
-------

Open and issue in the `issue tracker`_ for all support requests.

Reporting Issues
----------------

Report all issues in the `issue tracker`_. When doing so, please include
version information for:

- Python
- `py-smps`
- Operating System

Submitting Patches
------------------

All patches should be submitted as pull requests on the `GitHub project`_.

- Include tests if fixing a bug

- Clearly explain what you're trying to accomplish

- Follow :pep:`8`. You can use the `pep8` tox target for this

Testing
-------

`quantaq-cli` uses `coverage` and `unittest` for testing. To run all tests, run:

.. code-block:: shell

    $ poetry run coverage run -m unittest discover

We support a number of Python versions. Currently, we support and test for 3.6, 
3.7, and 3.8. These tests are automated and run as GitHub actions on every 
pull request.

.. _issue tracker: https://github.com/quant-aq/py-smps/issues
.. _GitHub project: https://github.com/quant-aq/py-smps