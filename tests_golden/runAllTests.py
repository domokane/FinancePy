###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("./..")

from os.path import dirname, basename, join
import glob
from financepy.utils.error import FinError
from financepy.utils.date import set_date_format, DateFormatTypes


# This only works if I have an init.py in the parent folder


set_date_format(DateFormatTypes.UK_LONG)

# I put this here to get the library loaded and header printed before loop

print("Looking in folder:", dirname(__file__))
modules = sorted(glob.glob(join(dirname(__file__), "Test*.py")))
numModules = len(modules)

""" This is the index of the file - change this to start later in the list """
n = 0
m = numModules

###############################################################################

for moduleFileName in modules[n:m+1]:

    try:

        moduleTextName = basename(moduleFileName[:-3])
        print("TEST: %3d out of %3d: MODULE: %-35s " % (n+1, numModules,
                                                        moduleTextName), end="")
        moduleName = __import__(moduleTextName)
        numErrors = moduleName.testCases._globalNumErrors
        numWarnings = moduleName.testCases._globalNumWarnings

        print("WARNINGS: %3d ERRORS: %3d " % (numWarnings, numErrors), end="")

        if numErrors > 0:
            for i in range(0, numErrors):
                print("*", end="")

        print("")
        n = n + 1

    # Want testing to continue even if a module has an exception
    except FinError as err:
        print("FinError:", err._message, "************")
        n = n + 1
        pass
    except ValueError as err:
        print("Value Error:", err.args[0], "************")
        n = n + 1
        pass
    except NameError as err:
        print("Name Error:", err.args[0], "************")
        n = n + 1
        pass
    except TypeError as err:
        print("Type Error:", err.args[0], "************")
        n = n + 1
        pass
    except BaseException as e:
        print("Base error:", e)
        n = n + 1
        pass
    except:
        print("Unexpected error:", sys.exc_info()[0])
        n = n + 1
        pass

###############################################################################
