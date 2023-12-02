###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys

sys.path.append("./..")

from financepy.utils.date import set_date_format, DateFormatTypes
from financepy.utils.error import FinError
import glob
from os.path import dirname, basename, join

# This only works if I have an init.py in the parent folder

set_date_format(DateFormatTypes.UK_LONG)

# I put this here to get the library loaded and header printed before loop

print("Looking in folder:", dirname(__file__))
modules = sorted(glob.glob(join(dirname(__file__), "Test*.py")))
num_modules = len(modules)

""" This is the index of the file - change this to start later in the list """
N = 0
M = num_modules

###############################################################################

for module_file_name in modules[N:M + 1]:

    try:

        module_text_name = basename(module_file_name[:-3])
        print("TEST: %3d out of %3d: MODULE: %-35s " % (N + 1, num_modules,
                                                        module_text_name),
              end="")
        module_name = __import__(module_text_name)
        num_errors = module_name.testCases._globalNumErrors
        num_warnings = module_name.testCases._globalNumWarnings

        print("WARNINGS: %3d ERRORS: %3d " % (num_warnings, num_errors),
              end="")

        if num_errors > 0:
            for i in range(0, num_errors):
                print("*", end="")

        print("")
        N = N + 1

    # Want testing to continue even if a module has an exception
    except FinError as err:
        print("FinError:", err._message, "************")
        N = N + 1
        pass
    except ValueError as err:
        print("Value Error:", err.args[0], "************")
        N = N + 1
        pass
    except NameError as err:
        print("Name Error:", err.args[0], "************")
        N = N + 1
        pass
    except TypeError as err:
        print("Type Error:", err.args[0], "************")
        N = N + 1
        pass
    except BaseException as e:
        print("Base error:", e)
        N = N + 1
        pass
    except Exception:
        print("Unexpected error:", sys.exc_info()[0])
        N = N + 1
        pass

###############################################################################
