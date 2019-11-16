import glob
from os.path import dirname, join

import sys
sys.path.append("..//..")

print("Looking in folder:", dirname(__file__))
modules = glob.glob(join(dirname(__file__), "TestFin*.py"))

numModules = len(modules)
n = 0
for moduleFileName in modules[0:]:

    n = n + 1
    ll = moduleFileName.find("\\TestFin")
    moduleTextName = moduleFileName[ll + 1:-3]

    print("==================================================================")
    print("TEST CASE ANALYSIS OF MODULE: ", moduleTextName)
    print("Module %d out of %d: " % (n, numModules))
    moduleName = __import__(moduleTextName)
