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

for moduleFileName in modules[n:m+1]:

    moduleTextName = basename(moduleFileName[:-3])

    print("==================================================================")
    print("Testing module %s (%d out of %d)"% (moduleTextName, n+1, numModules))
    moduleName = __import__(moduleTextName)
    n = n + 1

