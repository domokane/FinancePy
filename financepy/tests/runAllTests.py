import glob
from os.path import dirname, join

import sys
sys.path.append("..//..")

print("Looking in folder:", dirname(__file__))
modules = glob.glob(join(dirname(__file__), "TestFin*.py"))

for moduleFileName in modules:

    l = moduleFileName.find("\\TestFin")
    moduleTextName = moduleFileName[l + 1:-3]

    print("====================================================================")
    print("TEST CASE ANALYSIS OF MODULE: ", moduleTextName)
    moduleName = __import__(moduleTextName)
