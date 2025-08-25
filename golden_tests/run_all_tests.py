########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

import glob
from os.path import dirname, basename, join
import traceback

import sys
import os

# Add the *parent of the script's parent directory* to sys.path
sys.path.append(
    os.path.dirname(  # go up one directory
        os.path.dirname(  # go up two directories
            os.path.abspath(__file__)  # absolute path to current file
        )
    )
)

from financepy.utils.date_format import set_date_format, DateFormatTypes
from financepy.utils.error import FinError


# This only works if I have an init.py in the parent folder

set_date_format(DateFormatTypes.UK_LONG)


def main(start_index=0, end_index=None):
    """Loop over test cases"""

    # I put this here to get the library loaded and header printed before loop
    test_folder = dirname(__file__)
    print("Looking in folder:", test_folder)
    modules = sorted(glob.glob(join(test_folder, "Test*.py")))
    num_modules = len(modules)

    if end_index is None or end_index > num_modules:
        end_index = num_modules

    for idx in range(start_index, end_index):

        module_path = modules[idx]
        module_name = basename(module_path)[:-3]

        print(
            f"TEST: {idx + 1:3d} out of {num_modules:3d}: MODULE: {module_name:<35} ",
            end="",
        )

        try:
            module = __import__(module_name)

            num_errors = getattr(module.test_cases, "_global_num_errors", 0)
            num_warnings = getattr(module.test_cases, "_global_num_warnings", 0)

            print(f"WARNINGS: {num_warnings:3d} ERRORS: {num_errors:3d} ", end="")

            if num_errors > 0:
                print("*" * num_errors, end="")

            print()

        except (FinError, ValueError, NameError, TypeError) as e:
            # Handle known financepy or Python errors gracefully
            print(f"{type(e).__name__}: {e} ************")
        except Exception as e:
            # Catch all other unexpected errors, print traceback for debugging
            print(f"Unexpected {type(e).__name__}: {e}")
            traceback.print_exc()


if __name__ == "__main__":
    # Optionally customize start and end test indices here
    main()

########################################################################################
