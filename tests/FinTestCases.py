##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import sys
sys.path.append("..")

from os.path import join, exists, split
import time

from enum import Enum
from financepy.finutils.FinError import FinError

class FinTestCaseMode(Enum):
    SAVE_TEST_CASES = 1
    ANALYSE_TEST_CASES = 2
    DEBUG_TEST_CASES = 3


verbose = False

###############################################################################
###############################################################################

# THIS IS WHERE YOU CHANGE THE SETTING FOR RUNNING TEST CASES
# TO SAVE A GOLDEN FILE JUST COMMENT OUT SECOND LINE THEN UNCOMMENT IT FOR
# TESTING


globalTestCaseMode = FinTestCaseMode.SAVE_TEST_CASES
globalTestCaseMode = FinTestCaseMode.ANALYSE_TEST_CASES
#globalTestCaseMode = FinTestCaseMode.DEBUG_TEST_CASES

tolerance = 1e-8


###############################################################################
###############################################################################


class FinTestCases():
    ''' Test case framework for FinancePy.
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
    to change without triggering an error.'''

###############################################################################

    def __init__(self, moduleName, mode):
        ''' Create the TestCase given the module name and whether we are in
        GOLDEN or COMPARE mode. '''

        rootFolder, moduleFilename = split(moduleName)

        self._carefulMode = False
        self._verbose = False

        if mode in FinTestCaseMode:
            self._mode = mode
        else:
            raise FinError("Unknown TestCase Mode")

        if mode == FinTestCaseMode.DEBUG_TEST_CASES:
            # Don't do anything
            self._verbose = True
            return

        self._moduleName = moduleFilename[0:-3]
        self._foldersExist = True
        self._rootFolder = rootFolder
        self._headerFields = None
        self._globalNumWarnings = 0
        self._globalNumErrors = 0

#        print("Root folder:",self._rootFolder)
#        print("Modulename:",self._moduleName)

        self._goldenFolder = join(rootFolder, "golden")
        self._differencesFolder = join(rootFolder, "differences")

        if exists(self._goldenFolder) is False:
            print("Looking for:", self._goldenFolder)
            print("GOLDEN Folder DOES NOT EXIST. You must create it. Exiting")
            self._foldersExist = False
            return None

        self._compareFolder = join(rootFolder, "compare")

        if exists(self._compareFolder) is False:
            print("Looking for:", self._compareFolder)
            print("COMPARE Folder DOES NOT EXIST. You must create it. Exiting")
            self._foldersExist = False
            return None

        self._goldenFilename = join(self._goldenFolder,
                                    self._moduleName + "_GOLDEN.testLog")

        self._compareFilename = join(self._compareFolder,
                                     self._moduleName + "_COMPARE.testLog")

        self._differencesFilename = join(self._differencesFolder,
                                         self._moduleName + "_DIFFS.testLog")

        if self._mode == FinTestCaseMode.SAVE_TEST_CASES:

            print("GOLDEN Test Case Creation for module:", moduleFilename)

            if exists(self._goldenFilename) and self._carefulMode:
                overwrite = input("File " + self._goldenFilename
                                  + " exists. Overwrite (Y/N) ?")
                if overwrite == "N":
                    print("Not overwriting. Saving test cases failed.")
                    return
                elif overwrite == "Y":
                    print("Overwriting existing file...")

            creationTime = time.strftime("%Y%m%d_%H%M%S")
            f = open(self._goldenFilename, 'w')
            f.write("File Created on:" + creationTime)
            f.write("\n")
            f.close()

        else:

#            print("GENERATING NEW OUTPUT FOR MODULE", moduleFilename,
#                "FOR COMPARISON.")

            if exists(self._compareFilename) and self._carefulMode:
                overwrite = input("File " + self._compareFilename +
                                  " exists. Overwrite (Y/N) ?")
                if overwrite == "N":
                    print("Not overwriting. Saving test cases failed.")
                    return
                elif overwrite == "Y":
                    print("Overwriting existing file...")

#            print("Creating empty file",self._compareFilename)
            creationTime = time.strftime("%Y%m%d_%H%M%S")
            f = open(self._compareFilename, 'w')
            f.write("File Created on:" + creationTime)
            f.write("\n")
            f.close()

###############################################################################

    def print(self, *args):
        ''' Print comma separated output to GOLDEN or COMPARE directory. '''

        if self._mode == FinTestCaseMode.DEBUG_TEST_CASES:
            print(args)
            return

        if not self._foldersExist:
            print("Cannot print as GOLDEN and COMPARE folders don't exist")
            return

        if self._headerFields is None:
            print("ERROR: Need to set header fields before printing results")
        elif len(self._headerFields) != len(args):
            n1 = len(self._headerFields)
            n2 = len(args)
            raise FinError("ERROR: Number of data columns is " + str(n1)
                           + " but must equal " + str(n2)
                           + " to align with headers.")

        if self._mode == FinTestCaseMode.SAVE_TEST_CASES:
            filename = self._goldenFilename
        else:
            filename = self._compareFilename

        f = open(filename, 'a')
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
        ''' Print a banner on a line to the GOLDEN or COMPARE directory. '''

        if self._mode == FinTestCaseMode.DEBUG_TEST_CASES:
            print(txt)
            return

        if not self._foldersExist:
            print("Cannot print as GOLDEN and COMPARE folders do not exist")
            return

        if self._mode == FinTestCaseMode.SAVE_TEST_CASES:
            filename = self._goldenFilename
        else:
            filename = self._compareFilename

        f = open(filename, 'a')
        f.write("BANNER,")

        f.write(txt)

        f.write("\n")
        f.close()

###############################################################################

    def header(self, *args):
        ''' Print a header on a line to the GOLDEN or COMPARE directory. '''

        if self._mode == FinTestCaseMode.DEBUG_TEST_CASES:
            self.printLog(args)
            return

        if not self._foldersExist:
            self.printLog("Cannot print as GOLDEN and COMPARE folders do not exist")
            return

        self._headerFields = args

        if len(self._headerFields) == 0:
            self.printLog("ERROR: Number of header fields must be greater than 0")

        if self._mode == FinTestCaseMode.SAVE_TEST_CASES:
            filename = self._goldenFilename
        else:
            filename = self._compareFilename

        f = open(filename, 'a')
        f.write("HEADER,")

        for arg in args:
            f.write(str(arg))
            f.write(",")

        f.write("\n")
        f.close()

###############################################################################

    def compareRows(self, goldenRow, compareRow, rowNum):
        ''' Compare the contents of two rows in GOLDEN and COMPARE folders.'''

        numWarnings = 0
        numErrors = 0

        goldenFields = goldenRow.split(",")
        compareFields = compareRow.split(",")

        numGoldenFields = len(goldenFields)
        numCompareFields = len(compareFields)

        if numGoldenFields == 0:
            self.printLog("ERROR: No field in golden row")
            numErrors += 1
            return (0, numErrors)

        if numCompareFields == 0:
            self.printLog("ERROR: No field in golden row")
            numErrors += 1
            return (0, numErrors)

        if numGoldenFields != numCompareFields:
            self.printLog("ERROR: Mismatch in number of fields")
            numErrors += 1
            return (0, numErrors)

        numCols = numGoldenFields

        if goldenFields[0] != compareFields[0]:
            self.printLog("ERROR: Mismatch in row types HEADER vs RESULT",
                     goldenFields[0], compareFields[0])
            numErrors += 1
            return (0, numErrors)

        if goldenFields[0] == "HEADER" and compareFields[0] == "HEADER":
            for colNum in range(0, numCols):
                if compareFields[colNum] != goldenFields[colNum]:
                    numErrors += 1

            self._headerFields = goldenFields
            return (0, numErrors)

        if self._headerFields is None:
            self.printLog("ERROR: No Header row has been assigned")
            numErrors += 1
            return (0, numErrors)

        for colNum in range(0, numCols):

            if compareFields[colNum] != goldenFields[colNum]:

                if len(self._headerFields) <= colNum:
                    self.printLog("ERROR: Mismatch in headers. Rerun GOLDEN!")
                    return (0, numErrors)

                if self._headerFields[colNum] == "TIME":
                    time1 = float(goldenFields[colNum])
                    
#                    if isinstance(compareFields[colNum], str):
#                        self.printLog("Row# ", rowNum,
#                              " WARNING: The field in the new log is not a time")
                        

                    time2 = float(compareFields[colNum])
                    change = (time2 / abs(time1 + 1e-10) - 1.0) * 100.0

                    if abs(change) > 20.0:
                        self.printLog("Row# ", rowNum,
                              " WARNING: Calculation time has changed by %5.2f"
                              % change, " percent.")

                    numWarnings += 1
                else:
                    numErrors += 1

        return (numWarnings, numErrors)

###############################################################################

    def compareTestCases(self):
        ''' Compare output of COMPARE mode to GOLDEN output '''
        
        self.startLog()

        if self._mode == FinTestCaseMode.SAVE_TEST_CASES:
            # Do nothing
            return

        if self._mode == FinTestCaseMode.DEBUG_TEST_CASES:
            # Do nothing
            return

        self.printLog("EXAMINING CHANGES IN NEW OUTPUT AND GOLDEN FOR Module: " +
                 self._moduleName)

        totalNumWarnings = 0
        totalNumErrors = 0

        # check golden file exists
        if exists(self._goldenFilename) is False:
            print("No GOLDEN file exists. You need to change the mode to SAVE_TEST_CASES and then rerun.")
            return

        # open golden file and load it up
        f = open(self._goldenFilename, 'r')
        goldenContents = f.readlines()
        f.close()

        numGoldenLines = len(goldenContents)

        # open golden file and load it up
        f = open(self._compareFilename, 'r')
        compareContents = f.readlines()
        f.close()

        numCompareLines = len(compareContents)

        if numGoldenLines != numCompareLines:
            self.printLog("File lengths not the same")
            self.printLog("Number of COMPARE lines: ", numCompareLines)
            self.printLog("Number of GOLDEN lines: ", numGoldenLines)

        minNumLines = min(numGoldenLines, numCompareLines)
#        maxNumLines = max(numGoldenLines, numCompareLines)

        # We start at second row as first row has time stamp
        for rowNum in range(1, minNumLines):

            goldenRow = goldenContents[rowNum]
            compareRow = compareContents[rowNum]

            numWarnings, numErrors = self.compareRows(
                goldenRow, compareRow, rowNum)

            if numErrors > 0:
                self.printLog(
                    "Row# ",
                    rowNum,
                    " ERROR - ",
                    numErrors,
                    " DIFFERENCE(S) IN ROW ",
                    rowNum)

                if self._headerFields is None:
                    self.printLog("ERROR: Header must be defined")
                    numErrors += 1
                else:
                    self.printLog("Row# ", rowNum, " HEADER:  ==>",
                             self._headerFields[0:-1])

                self.printLog("Row# ", rowNum, " GOLDEN : ==>", goldenRow[:-2])
                self.printLog("Row# ", rowNum, " COMPARE: ==>", compareRow[:-2])
                self.printLog("")

            totalNumWarnings += numWarnings
            totalNumErrors += numErrors

        if numGoldenLines == minNumLines and numCompareLines > minNumLines:
            self.printLog("ERROR:The COMPARE file is longer than the GOLDEN file")

#            for rowNum in range(minNumLines,maxNumLines):
#                numWarnings, numErrors = compareContents[rowNum]
#                print(rowNum,"COMPARE: ==>",compareRow)

        if numCompareLines == minNumLines and numGoldenLines > minNumLines:
            self.printLog("ERROR:The GOLDEN file is longer than the COMPARE file")

#            for rowNum in range(minNumLines,maxNumLines):
#                numWarnings, numErrors = compareContents[rowNum]
#                print(rowNum,"GOLDEN: ==>",goldenRow)

#        print("Analysis of", self._moduleName, "completed with",
#              totalNumErrors, "errors and", totalNumWarnings, "warnings.")

#        print("NUM LINES:", numCompareLines, 
#              "====>",
#              "ERRORS:", totalNumErrors, 
#              "WARNINGS:", totalNumWarnings)

        self.printLog("Analysis of ", self._moduleName, " completed with ",
                      totalNumErrors, " errors and ", totalNumWarnings,
                      " warnings.")

        self.printLog("NUM LINES:", numCompareLines,
                      "====>",
                      "ERRORS:", totalNumErrors, 
                      " WARNINGS:", totalNumWarnings)

        self._globalNumErrors = totalNumErrors
        self._globalNumWarnings = totalNumWarnings

        return

###############################################################################

    def startLog(self):
        f = open(self._differencesFilename, 'w')
        f.close()

###############################################################################

    def printLog(self, *args):
        f = open(self._differencesFilename, 'a')
        for arg in args:
            f.write(str(arg) + " ")

        f.write("\n")
        f.close()

###############################################################################
###############################################################################
