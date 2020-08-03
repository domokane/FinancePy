This is a collection of modules used across a wide range of FinancePy functions. Examples include date generation, special mathematical functions and useful helper functions for performing some repeated action.

* FinDate is a class for handling dates in a financial setting. Special functions are included for computing IMM dates and CDS dates and moving dates forward by tenors.
* FinCalendar is a class for determining which dates are not business dates in a specific region or country.
* FinDayCount is a class for determining accrued interest in bonds and also accrual factors in ISDA swap-like contracts.
* FinError is a class which handles errors in the calculations done within FinancePy
* FinFrequency takes in a frequency type and then returns the number of payments per year
* FinGlobalVariables holds the value of constants used across the whole of FinancePy
* FinHelperFunctions is a set of helpful functions that can be used in a number of places
* FinMath is a set of mathematical functions specific to finance which have been optimised for speed using Numba
* FinSobol is the implementation of Sobol quasi-random number generator. It has been speeded up using Numba.
* FinRateConverter converts rates for one compounding frequency to rates for a different frequency
* FinSchedule generates a sequence of cashflow payment dates in accordance with financial market standards
* FinStatistics calculates a number of statistical variables such as mean, standard deviation and variance
* FinTestCases is the code that underlies the test case framework used across FinancePy
