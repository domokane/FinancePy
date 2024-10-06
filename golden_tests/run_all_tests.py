###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import glob
from os.path import dirname, basename, join

import sys

sys.path.append("..")

from financepy.utils.date import set_date_format, DateFormatTypes
from financepy.utils.error import FinError

# This only works if I have an init.py in the parent folder

set_date_format(DateFormatTypes.UK_LONG)

# I put this here to get the library loaded and header printed before loop

print("Looking in folder:", dirname(__file__))
modules = sorted(glob.glob(join(dirname(__file__), "Test*.py")))
num_modules = len(modules)

""" This is the index of the file - change this to start later in the list """
start_module = 0
end_module = num_modules

###############################################################################

module_number = start_module

for module_file_name in modules[start_module : end_module + 1]:

    try:

        module_text_name = basename(module_file_name[:-3])
        print(
            "TEST: %3d out of %3d: MODULE: %-35s "
            % (module_number + 1, num_modules, module_text_name),
            end="",
        )

        module_name = __import__(module_text_name)
        num_errors = module_name.test_cases._global_num_errors
        num_warnings = module_name.test_cases._global_num_warnings

        print(
            "WARNINGS: %3d ERRORS: %3d " % (num_warnings, num_errors), end=""
        )

        if num_errors > 0:
            for i in range(0, num_errors):
                print("*", end="")

        print("")
        module_number = module_number + 1

    # Want testing to continue even if a module has an exception
    except FinError as err:
        print("FinError:", err._message, "************")
        module_number = module_number + 1
    except ValueError as err:
        print("Value Error:", err.args[0], "************")
        module_number = module_number + 1
    except NameError as err:
        print("Name Error:", err.args[0], "************")
        module_number = module_number + 1
    except TypeError as err:
        print("Type Error:", err.args[0], "************")
        module_number = module_number + 1
    except Exception:
        print("Unexpected error:", sys.exc_info()[0])
        module_number = module_number + 1

###############################################################################
