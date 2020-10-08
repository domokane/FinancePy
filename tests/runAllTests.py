###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import glob
from os.path import dirname, basename, join

import sys
sys.path.append("..")

print("Looking in folder:", dirname(__file__))
modules = sorted(glob.glob(join(dirname(__file__), "Test*.py")))

numModules = len(modules)

''' This is the index of the file - change this to start later in the list '''
n = 0
m = numModules

# I put this here to get the header print before the loop starts
from FinTestCases import FinTestCases

for moduleFileName in modules[n:m+1]:

    moduleTextName = basename(moduleFileName[:-3])

    print("TEST: %3d out of %3d: MODULE: %-35s "% (n+1, numModules, moduleTextName), end="")
    moduleName = __import__(moduleTextName)

    numErrors = moduleName.testCases._globalNumErrors
    numWarnings = moduleName.testCases._globalNumWarnings

    print("WARNINGS: %3d   ERRORS: %3d " % (numWarnings, numErrors), end ="")

    if numErrors > 0:
        for i in range(0, numErrors):
            print("*", end="")
    
    print("")

    n = n + 1

