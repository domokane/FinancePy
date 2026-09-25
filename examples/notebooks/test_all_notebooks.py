# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 12:36:01 2021

@author: Dominic
"""

import asyncio
import glob
import os
import sys

import nbformat
from nbconvert.preprocessors import CellExecutionError, ExecutePreprocessor

if sys.platform.startswith("win"):
    asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())


print("Looking in folder:", os.path.dirname(__file__))

notebooks = sorted(
    glob.glob(
        os.path.join(os.path.dirname(__file__), "./products/*/*.ipynb"), recursive=True
    )
)


########################################################################################


def notebook_run_new(notebook_filepathname: str) -> None:
    """Execute a notebook via nbconvert and collect output."""

    notebook_dir, filename = os.path.split(notebook_filepathname)

    os.chdir(notebook_dir)

    with open(filename, encoding="utf-8") as file:
        nb = nbformat.read(file, as_version=4)

    ep = ExecutePreprocessor(timeout=600, kernel_name="python3")

    notebook_filename_err = os.path.join(notebook_dir, "ERROR_" + filename)

    try:
        ep.preprocess(nb, {"metadata": {"path": "./"}})

        # Save notebook
        with open(notebook_filepathname, mode="w", encoding="utf-8") as file:
            nbformat.write(nb, file)

    except CellExecutionError:
        msg = f'Error executing the notebook "{filename}".\n'
        print(msg)

        with open(notebook_filename_err, mode="w", encoding="utf-8") as file:
            nbformat.write(nb, file)


########################################################################################


print("Starting")

n = 0
m = len(notebooks)

print(n, m)

for notebook in notebooks[n : m + 1]:
    notebook_dir, filename = os.path.split(notebook)

    print("Checking Notebook", n + 1, "of", m, ":", filename)

    notebook_run_new(notebook)
    n += 1
