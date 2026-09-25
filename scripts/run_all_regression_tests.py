########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

from financepy.utils.error import FinError
from financepy.utils.date_format import DateFormatTypes, set_date_format
import financepy
import glob
import sys
import time
import traceback
from os.path import basename, join
from pathlib import Path

# Ensure we import FinancePy from this repository
repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root))

print("These tests run against local repo and not installed FinancePy. They are mainly used for local development work.\n")

print("FinancePy imported successfully from:", financepy.__file__)
print()

set_date_format(DateFormatTypes.UK_LONG)

###############################################################################


def main(start_index=0, end_index=None):
    scripts_dir = Path(__file__).resolve().parent
    test_folder = scripts_dir.parent / "tests" / "regression"

    sys.path.insert(0, str(test_folder))

    modules = sorted(glob.glob(join(test_folder, "Test*.py")))

    num_modules = len(modules)

    if end_index is None or end_index > num_modules:
        end_index = num_modules

    timings = []

    # start_index = 40
    #    end_index = 115

    for idx in range(start_index, end_index):

        module_path = modules[idx]
        module_name = basename(module_path)[:-3]

        print(
            f"TEST: {idx + 1:3d} of {num_modules:3d}: {module_name:<35} ",
            end="",
        )

        start_time = time.perf_counter()

        try:
            module = __import__(module_name)

            num_errors = getattr(module.test_cases, "_global_num_errors", 0)
            num_warnings = getattr(module.test_cases, "_global_num_warnings", 0)

            elapsed = time.perf_counter() - start_time
            timings.append((module_name, elapsed))

            # print(f"WARNINGS: {num_warnings:3d} ERRORS: {num_errors:3d} ", end="")
            print(
                f"TIME: {elapsed:6.3f} s " f"WARNINGS: {num_warnings:3d} ERRORS: {num_errors:3d}",
                end="",
            )

            if num_errors > 0:
                print("*" * 1, end="")

            print()

        except (FinError, ValueError, NameError, TypeError) as e:
            elapsed = time.perf_counter() - start_time
            timings.append((module_name, elapsed))
            print(f"{type(e).__name__}: {e} ************ (TIME: {elapsed:6.3f} s)")
        except Exception as e:
            elapsed = time.perf_counter() - start_time
            timings.append((module_name, elapsed))
            print(f"Unexpected {type(e).__name__}: {e} (TIME: {elapsed:6.3f} s)")
            traceback.print_exc()


if __name__ == "__main__":
    # Optionally customize start and end test indices here
    main()

###############################################################################
