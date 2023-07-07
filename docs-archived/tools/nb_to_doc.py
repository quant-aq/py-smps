#! /usr/bin/env python
"""
Convert empty IPython notebook to a sphinx doc page.
"""
import sys
from subprocess import check_call as sh


def convert_nb(nbname):

    # Execute the notebook
    sh(["jupyter", "nbconvert", "--to", "notebook",
        "--execute", "--inplace", nbname + ".ipynb", 
        "--ExecutePreprocessor.timeout=60"])

    # Convert to .rst for Sphinx
    sh(["jupyter", "nbconvert", "--to", "rst", nbname + ".ipynb",
        "--ExecutePreprocessor.timeout=60"])

    # Clear notebook output
    sh(["jupyter", "nbconvert", "--to", "notebook", "--inplace",
        "--ClearOutputPreprocessor.enabled=True", nbname + ".ipynb",
        "--ExecutePreprocessor.timeout=60"])


if __name__ == "__main__":

    for nbname in sys.argv[1:]:
        convert_nb(nbname)
