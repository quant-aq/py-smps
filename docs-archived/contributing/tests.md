# Tests

_py-smps_ uses _coverage_ and _unittest_ for testing. To run all tests:

```sh
$ poetry run coverage run -m unittest discover
```

We support a number of python versions (Python 3.7+). Currently, we test against Python3.7, Python3.8, and Python3.9. These tests are automated and run via GitHub actions on each pull request.