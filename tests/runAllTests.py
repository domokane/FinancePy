import glob
from os.path import dirname

import sys
sys.path.append("..//..")

print("Looking in folder:", dirname(__file__))
modules = glob.glob("TestFin*.py")

numModules = len(modules)

''' This is the index of the file - change this to start later in the list '''
n = 69


for moduleFileName in modules[n:]:

    n = n + 1
    ll = moduleFileName.find("\\TestFin")
    moduleTextName = moduleFileName[ll + 1:-3]

    print("==================================================================")
    print("TEST CASE ANALYSIS OF MODULE: ", moduleTextName)
    print("Module %d out of %d: " % (n, numModules))
    moduleName = __import__(moduleTextName)
