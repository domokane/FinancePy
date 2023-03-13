# Test Suite

This folder contains the logic for performing comparison testing of code to ensure that code does not inadvertently get broken. 

The code tests here can be very broad in type. 

The test is based on the prior creation of a Golden version of the file output which is stored in the Golden folder. There is a file for each test file. This is considered to be correct and any output that deviates from this is deemed to be incorrect. Creation of the Golden files is done by setting the value of the global flag as follows:

globalTestCaseMode = FinTestCaseMode.SAVE_TEST_CASES

This can be found under utils/FinTestCase.py

When you wish to run the tests to check if anything has changed, you need to set the global flag to

globalTestCaseMode = FinTestCaseMode.ANALYSE_TEST_CASES

## RunAllTests.py

This file executes all of the test files which each report on any errors found.


## Output

Any change is registered as an error. No matter what the precision. 