###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import glob
from os.path import dirname, basename, join

import sys
sys.path.append("..")

print("Looking in folder:", dirname(__file__))
modules = glob.glob(join(dirname(__file__), "TestFin*.py"))

numModules = len(modules)

''' This is the index of the file - change this to start later in the list '''
n = 0
m = -1

for moduleFileName in modules[n:m]:

    n = n + 1
    moduleTextName = basename(moduleFileName[:-3])

    print("==================================================================")
    print("TEST CASE ANALYSIS OF MODULE: ", moduleTextName)
    print("Module %d out of %d: " % (n, numModules))
    moduleName = __import__(moduleTextName)
