# FinancePy Golden Regression Test Suite

This folder contains the logic for performing comparison testing of model-related code.

Because model-related code may produce slightly different results if a library dependency changes, I have not used the pytest framework here.

Also, I wish to test all of the parameter dimensions to check the full range of inputs.

I therefore use my own golden tests environment. This gathers the output and then compares it to a GOLDEN output set to examine changes.

You can examine any differences to see if they are problematic or not.

If all is OK. Just make the current output the new GOLDEN output.

## How to Code the Test Cases
Insert the following at the top of the test case file

    from FinTestCases import FinTestCases, global_test_case_mode

    test_cases = FinTestCases(__file__, global_test_case_mode)


Then in the last line of the module add the line

    test_cases.compare_test_cases()

To create values, you first define a header for the output.

### Example
This example header has 10 columns of model output that the tests will produce.

    test_cases.header(
       "STRIKE",
       "STEPS",
       "CALL_INT",
       "CALL_INT_PV",
       "CALL_EUR",
       "CALL_AMER",
       "PUT_INT",
       "PUT_INT_PV",
       "PUT_EUR",
       "PUT_AMER")
    ....
    .... calling model code to get values for all of these variables
    ....
    test_cases.print(strike_price,
                     num_steps,
                     call_intrinsic,
                     call_intrinsic_pv,
                     v1,
                     v2,
                     put_intrinsic,
                     put_intrinsic_pv,
                     v3,
                     v4,)

When you execute ***run_all_tests.py***, each test is run in order and a one-line summary of the results are output for each test file.

You can see WARNINGS which reflect timing changes and ERRORS which signify value changes

To investigate these warnings and errors, look at the output that goes to the data files that will be used to examine model output.

- To see the GOLDEN model output look in the *regression_tests/golden* subfolder here
- To see the latest model output look in the *regression_tests/compare* subfolder here
- To see the differencesm, look in the *regression_tests/differences* subfolder here

Each sub-folder contains a file for each test file that you can examine.

***THESE FILES ARE ACCESSING THE CURRENT CODE IN THE FINANCEPY FOLDER SO REGRESSION TESTING CAN BE USED AS PART OF THE MODEL DEVELOPMENT AND TESTING PROCESS***

### Creating the Golden files

The test is based on the prior creation of a Golden version of the file output which is stored in the Golden folder. There is a file for each test file. This is considered to be correct and any output that deviates from this is deemed to be incorrect. Creation of the Golden files is done by setting the value of the global flag as follows:

    globalTestCaseMode = FinTestCaseMode.SAVE_TEST_CASES

This can be found under FinTestCase.py

When you wish to run the tests to check if anything has changed, you need to set the global flag to

    globalTestCaseMode = FinTestCaseMode.ANALYSE_TEST_CASES

## Output

The framework compares the latest generated output against the corresponding GOLDEN file.

Numerical differences larger than the configured tolerance are reported as ERRORS. Timing columns labelled `TIME` are treated specially and may generate WARNINGS rather than ERRORS.

Warnings and errors should be inspected manually. If the new output is judged to be correct, the GOLDEN files can be regenerated.
