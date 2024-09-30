##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import sys

sys.path.append("..")

from financepy.utils.error import FinError
from enum import Enum
import time
from os.path import join, exists, split


class FinTestCaseMode(Enum):
    SAVE_TEST_CASES = 1
    ANALYSE_TEST_CASES = 2
    DEBUG_TEST_CASES = 3


VERBOSE = False

###############################################################################
###############################################################################

# THIS IS WHERE YOU CHANGE THE SETTING FOR RUNNING TEST CASES
# TO SAVE A GOLDEN FILE JUST COMMENT OUT SECOND LINE THEN UNCOMMENT IT FOR
# TESTING


# globalTestCaseMode = FinTestCaseMode.SAVE_TEST_CASES
globalTestCaseMode = FinTestCaseMode.ANALYSE_TEST_CASES
# globalTestCaseMode = FinTestCaseMode.DEBUG_TEST_CASES

TOLERANCE = 1e-8


###############################################################################
###############################################################################


class FinTestCases:
    """Test case framework for FinancePy.
    - The basic step is that we generate a GOLDEN folder that creates an output
    file for each testcase which is assumed to be correct. This can be done by
    running the test cases Python file with the globalTestCaseMode flag set to
    FinTestCaseMode.SAVE_TEST_CASES.
    - The second step is that we change the value of globalTestCaseMode to
    FinTestCaseMode.ANALYSE_TEST_CASES and then run the test scripts. This time
    they save a copy of the output to the COMPARE folder.
    Finally, a function called compareTestCases() is used to compare the new
    output with the GOLDEN output and states whether anything has changed.

    - The output of a test case has three forms each with its own method:

        1) print - this outputs comma separated values
        2) header - this must precede any print statement and labels the output
           columns
        3) banner - this is any single string line separator

    Note that the header TIME is special as it tells the analysis that the
    value in the corresponding column is a timing and so its value is allowed
    to change without triggering an error."""

    ###############################################################################

    def __init__(self, module_name, mode):
        """Create the TestCase given the module name and whether we are in
        GOLDEN or COMPARE mode."""

        root_folder, module_file_name = split(module_name)

        self._careful_mode = False
        self._verbose = False

        if mode in FinTestCaseMode:
            self._mode = mode
        else:
            raise FinError("Unknown TestCase Mode")

        if mode == FinTestCaseMode.DEBUG_TEST_CASES:
            # Don't do anything
            self._verbose = True
            return

        self._module_name = module_file_name[0:-3]
        self._foldersExist = True
        self._root_folder = root_folder
        self._header_fields = None
        self._global_num_warnings = 0
        self._global_num_errors = 0

        #        print("Root folder:",self._root_folder)
        #        print("module_name:",self._module_name)

        self._goldenFolder = join(root_folder, "golden")
        self._differencesFolder = join(root_folder, "differences")

        if exists(self._goldenFolder) is False:
            print("Looking for:", self._goldenFolder)
            print("GOLDEN Folder DOES NOT EXIST. You must create it. Exiting")
            self._foldersExist = False
            return

        self._compare_folder = join(root_folder, "compare")

        if exists(self._compare_folder) is False:
            print("Looking for:", self._compare_folder)
            print("COMPARE Folder DOES NOT EXIST. You must create it. Exiting")
            self._foldersExist = False
            return

        self._golden_file_name = join(
            self._goldenFolder, self._module_name + "_GOLDEN.testLog"
        )

        self._compare_file_name = join(
            self._compare_folder, self._module_name + "_COMPARE.testLog"
        )

        self._differences_file_name = join(
            self._differencesFolder, self._module_name + "_DIFFS.testLog"
        )

        if self._mode == FinTestCaseMode.SAVE_TEST_CASES:

            print("GOLDEN Test Case Creation for module:", module_file_name)

            if exists(self._golden_file_name) and self._careful_mode:
                overwrite = input(
                    "File "
                    + self._golden_file_name
                    + " exists. Overwrite (Y/N) ?"
                )
                if overwrite == "N":
                    print("Not overwriting. Saving test cases failed.")
                    return
                elif overwrite == "Y":
                    print("Overwriting existing file...")

            creation_time = time.strftime("%Y%m%d_%H%M%S")
            f = open(self._golden_file_name, "w", encoding="utf-8")
            f.write("File Created on:" + creation_time)
            f.write("\n")
            f.close()

        else:

            #            print("GENERATING NEW OUTPUT FOR MODULE", module_file_name,
            #                "FOR COMPARISON.")

            if exists(self._compare_file_name) and self._careful_mode:
                overwrite = input(
                    "File "
                    + self._compare_file_name
                    + " exists. Overwrite (Y/N) ?"
                )
                if overwrite == "N":
                    print("Not overwriting. Saving test cases failed.")
                    return
                elif overwrite == "Y":
                    print("Overwriting existing file...")

            #            print("Creating empty file",self._compare_file_name)
            creation_time = time.strftime("%Y%m%d_%H%M%S")
            f = open(self._compare_file_name, "w", encoding="utf-8")
            f.write("File Created on:" + creation_time)
            f.write("\n")
            f.close()

    ###############################################################################

    def print(self, *args):
        """Print comma separated output to GOLDEN or COMPARE directory."""

        if self._mode == FinTestCaseMode.DEBUG_TEST_CASES:
            print(args)
            return

        if not self._foldersExist:
            print("Cannot print as GOLDEN and COMPARE folders don't exist")
            return

        if self._header_fields is None:
            print("ERROR: Need to set header fields before printing results")
        elif len(self._header_fields) != len(args):
            n1 = len(self._header_fields)
            n2 = len(args)
            raise FinError(
                "ERROR: Number of data columns is "
                + str(n1)
                + " but must equal "
                + str(n2)
                + " to align with headers."
            )

        if self._mode == FinTestCaseMode.SAVE_TEST_CASES:
            file_name = self._golden_file_name
        else:
            file_name = self._compare_file_name

        f = open(file_name, "a", encoding="utf-8")
        f.write("RESULTS,")

        for arg in args:
            if isinstance(arg, float):
                f.write("%10.8f" % (arg))
            else:
                f.write(str(arg))
            f.write(",")

        f.write("\n")
        f.close()

    ###############################################################################

    def banner(self, txt):
        """Print a banner on a line to the GOLDEN or COMPARE directory."""

        if self._mode == FinTestCaseMode.DEBUG_TEST_CASES:
            print(txt)
            return

        if not self._foldersExist:
            print("Cannot print as GOLDEN and COMPARE folders do not exist")
            return

        if self._mode == FinTestCaseMode.SAVE_TEST_CASES:
            file_name = self._golden_file_name
        else:
            file_name = self._compare_file_name

        f = open(file_name, "a", encoding="utf-8")
        f.write("BANNER,")

        f.write(txt)

        f.write("\n")
        f.close()

    ###############################################################################

    def header(self, *args):
        """Print a header on a line to the GOLDEN or COMPARE directory."""

        if self._mode == FinTestCaseMode.DEBUG_TEST_CASES:
            self.print_log(args)
            return

        if not self._foldersExist:
            self.print_log(
                "Cannot print as GOLDEN and COMPARE folders do not exist"
            )
            return

        self._header_fields = args

        if len(self._header_fields) == 0:
            self.print_log(
                "ERROR: Number of header fields must be greater than 0"
            )

        if self._mode == FinTestCaseMode.SAVE_TEST_CASES:
            file_name = self._golden_file_name
        else:
            file_name = self._compare_file_name

        f = open(file_name, "a", encoding="utf-8")
        f.write("HEADER,")

        for arg in args:
            f.write(str(arg))
            f.write(",")

        f.write("\n")
        f.close()

    ###############################################################################

    def compare_rows(self, golden_row, compare_row, row_num):
        """Compare the contents of two rows in GOLDEN and COMPARE folders."""

        num_warnings = 0
        num_errors = 0

        golden_fields = golden_row.split(",")
        compare_fields = compare_row.split(",")

        num_golden_fields = len(golden_fields)
        num_compare_fields = len(compare_fields)

        if num_golden_fields == 0:
            self.print_log("ERROR: No field in golden row")
            num_errors += 1
            return (0, num_errors)

        if num_compare_fields == 0:
            self.print_log("ERROR: No field in golden row")
            num_errors += 1
            return (0, num_errors)

        if num_golden_fields != num_compare_fields:
            self.print_log("ERROR: Mismatch in number of fields")
            num_errors += 1
            return (0, num_errors)

        num_cols = num_golden_fields

        if golden_fields[0] != compare_fields[0]:
            self.print_log(
                "ERROR: Mismatch in row types HEADER vs RESULT",
                golden_fields[0],
                compare_fields[0],
            )
            num_errors += 1
            return (0, num_errors)

        if golden_fields[0] == "HEADER" and compare_fields[0] == "HEADER":
            for col_num in range(0, num_cols):
                if compare_fields[col_num] != golden_fields[col_num]:
                    num_errors += 1

            self._header_fields = golden_fields
            return (0, num_errors)

        if self._header_fields is None:
            self.print_log("ERROR: No Header row has been assigned")
            num_errors += 1
            return (0, num_errors)

        for col_num in range(0, num_cols):

            if len(self._header_fields) <= col_num:
                self.print_log("ERROR: Mismatch in headers. Rerun!")
                return (0, num_errors)

            time_column = False
            if self._header_fields[col_num] == "TIME":
                time_column = True

            compare_field = compare_fields[col_num]
            golden_field = golden_fields[col_num]

            if compare_field != golden_field:

                # Only reject with TOLERANCE of 1e-8 - OK for finance!
                tol = 1e-4
                compare_flag = compare_field.replace(".", "").isnumeric()
                golden_flag = golden_field.replace(".", "").isnumeric()

                if (
                    compare_flag is True
                    and golden_flag is True
                    and time_column is False
                ):

                    compare_value = float(compare_field)
                    golden_value = float(golden_field)
                    err = compare_value - golden_value

                    if abs(err) > tol:
                        num_errors += 1
                #                        print("OK:", compare_value, golden_value, num_errors)

                if time_column is True:
                    time1 = float(golden_fields[col_num])
                    time2 = float(compare_fields[col_num])
                    change = (time2 / abs(time1 + 1e-10) - 1.0) * 100.0

                    if abs(change) > 50.0:
                        self.print_log(
                            "Row# ",
                            row_num,
                            " WARNING: Calculation time has changed by %5.2f"
                            % change,
                            " percent.",
                        )

                    num_warnings += 1

        return (num_warnings, num_errors)

    ###############################################################################

    def compareTestCases(self):
        """Compare output of COMPARE mode to GOLDEN output"""

        self.start_log()

        if self._mode == FinTestCaseMode.SAVE_TEST_CASES:
            # Do nothing
            return

        if self._mode == FinTestCaseMode.DEBUG_TEST_CASES:
            # Do nothing
            return

        self.print_log(
            "EXAMINING CHANGES IN NEW OUTPUT AND GOLDEN FOR Module: "
            + self._module_name
        )

        totalnum_warnings = 0
        totalnum_errors = 0

        # check golden file exists
        if exists(self._golden_file_name) is False:
            print(
                "No GOLDEN file exists. You need to change the mode to SAVE_TEST_CASES and then rerun."
            )
            return

        # open golden file and load it up
        with open(self._golden_file_name, "r", encoding="utf-8") as f:
            golden_contents = f.readlines()

        num_golden_lines = len(golden_contents)

        # open golden file and load it up
        with open(self._compare_file_name, "r", encoding="utf-8") as f:
            compareContents = f.readlines()

        num_compare_lines = len(compareContents)

        if num_golden_lines != num_compare_lines:
            self.print_log("File lengths not the same")
            self.print_log("Number of COMPARE lines: ", num_compare_lines)
            self.print_log("Number of GOLDEN lines: ", num_golden_lines)

        minNumLines = min(num_golden_lines, num_compare_lines)
        #        maxNumLines = max(num_golden_lines, num_compare_lines)

        # We start at second row as first row has time stamp
        for row_num in range(1, minNumLines):

            golden_row = golden_contents[row_num]
            compare_row = compareContents[row_num]

            num_warnings, num_errors = self.compare_rows(
                golden_row, compare_row, row_num
            )

            if num_errors > 0:
                self.print_log(
                    "Row# ",
                    row_num,
                    " ERROR - ",
                    num_errors,
                    " DIFFERENCE(S) IN ROW ",
                    row_num,
                )

                if self._header_fields is None:
                    self.print_log("ERROR: Header must be defined")
                    num_errors += 1
                else:
                    self.print_log(
                        "Row# ",
                        row_num,
                        " HEADER:  ==>",
                        self._header_fields[0:-1],
                    )

                self.print_log(
                    "Row# ", row_num, " GOLDEN : ==>", golden_row[:-2]
                )
                self.print_log(
                    "Row# ", row_num, " COMPARE: ==>", compare_row[:-2]
                )
                self.print_log("")

            totalnum_warnings += num_warnings
            totalnum_errors += num_errors

        if num_golden_lines == minNumLines and num_compare_lines > minNumLines:
            self.print_log(
                "ERROR:The COMPARE file is longer than the GOLDEN file"
            )

        #            for row_num in range(minNumLines,maxNumLines):
        #                num_warnings, num_errors = compareContents[row_num]
        #                print(row_num,"COMPARE: ==>",compare_row)

        if num_compare_lines == minNumLines and num_golden_lines > minNumLines:
            self.print_log(
                "ERROR:The GOLDEN file is longer than the COMPARE file"
            )

        #            for row_num in range(minNumLines,maxNumLines):
        #                num_warnings, num_errors = compareContents[row_num]
        #                print(row_num,"GOLDEN: ==>",golden_row)

        #        print("Analysis of", self._module_name, "completed with",
        #              totalnum_errors, "errors and", totalnum_warnings, "warnings.")

        #        print("NUM LINES:", num_compare_lines,
        #              "====>",
        #              "ERRORS:", totalnum_errors,
        #              "WARNINGS:", totalnum_warnings)

        self.print_log(
            "Analysis of ",
            self._module_name,
            " completed with ",
            totalnum_errors,
            " errors and ",
            totalnum_warnings,
            " warnings.",
        )

        self.print_log(
            "NUM LINES:",
            num_compare_lines,
            "====>",
            "ERRORS:",
            totalnum_errors,
            " WARNINGS:",
            totalnum_warnings,
        )

        self._global_num_errors = totalnum_errors
        self._global_num_warnings = totalnum_warnings

        return

    ###############################################################################

    def start_log(self):
        f = open(self._differences_file_name, "w", encoding="utf-8")
        f.close()

    ###############################################################################

    def print_log(self, *args):
        f = open(self._differences_file_name, "a", encoding="utf-8")
        for arg in args:
            f.write(str(arg) + " ")

        f.write("\n")
        f.close()


###############################################################################
###############################################################################
