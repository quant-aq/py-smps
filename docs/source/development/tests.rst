Tests
=====

`py-smps` uses the *coverage* and *unittest* libraries to run automated tests. To manually 
run unittests:

.. code-block:: console

    $ poetry run coverage run -m --source=. pytest tests/ --ignore=tests/datafiles

To view the coverage report, you can then run:

.. code-block:: console

    $ poetry run coverage report -m