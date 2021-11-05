# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 12:36:01 2021

@author: Dominic
"""

import os

import nbformat

from nbconvert.preprocessors import ExecutePreprocessor
from nbconvert.preprocessors import CellExecutionError

import glob
from os.path import dirname, join

print("Looking in folder:", dirname(__file__))

notebooks = sorted(
    glob.glob(join(dirname(__file__), "./products/*/*.ipynb"), recursive=True))

###############################################################################


def notebook_run_new(notebook_filepathname):
    """Execute a notebook via nbconvert and collect output.
       :returns (parsed nb object, execution errors)
    """

    dirname, filename = os.path.split(notebook_filepathname)

    os.chdir(dirname)

    nb = nbformat.read(open(notebook), as_version=4)
    ep = ExecutePreprocessor(timeout=600, kernel_name='python3')

    notebook_filename_out = dirname + "\ERROR_" + filename

    try:
        out = ep.preprocess(nb, {'metadata': {'path': "./"}})
    except CellExecutionError:

        out = None
        msg = 'Error executing the notebook "%s".\n\n' % filename
        msg += 'See notebook "%s" for the traceback.' % notebook_filename_out
        print(msg)

        with open(notebook_filename_out, mode='w', encoding='utf-8') as f:
            nbformat.write(nb, f)

        pass

###############################################################################


print("Starting")
n = 0
m = len(notebooks)
print(n, m)
for notebook in notebooks[n:m+1]:
    dirname, filename = os.path.split(notebook)
    print("Checking Notebook", filename, "which is", n+1, "of", m)
    notebook_run_new(notebook)
    n = n + 1
