# Introduction 

This is a collection of modules used across a wide range of FinancePy functions. Examples include date generation, special mathematical functions and useful helper functions for performing some repeated action

* Date is a class for handling dates in a financial setting. Special functions are included for computing IMM dates and CDS dates and moving dates forward by tenors.
* Calendar is a class for determining which dates are not business dates in a specific region or country.
* DayCount is a class for determining accrued interest in bonds and also accrual factors in ISDA swap-like contracts.
* Error is a class which handles errors in the calculations done within FinancePy
* annual_frequency takes in a annual_frequency type and then returns the number of payments per year
* global_vars holds the value of constants used across the whole of FinancePy
* helper_functions is a set of helpful functions that can be used in a number of places
* math is a set of mathematical functions specific to finance which have been optimised for speed using Numba
* FinSobol is the implementation of Sobol quasi-random number generator. It has been speeded up using Numba.
* FinRateConverter converts rates for one compounding annual_frequency to rates for a different annual_frequency
* FinSchedule generates a sequence of cashflow payment dates in accordance with financial market standards
* FinStatistics calculates a number of statistical variables such as mean, standard deviation and variance
* FinTestCases is the code that underlies the test case framework used across FinancePy

## FinDayCount

The year fraction function can take up to 3 dates, D1, D2 and D3 and a annual_frequency in specific cases. The current day count methods are listed below.

* THIRTY 360 BOND - 30E/360 ISDA 2006 4.16f, German, Eurobond(ISDA 2000)
* THIRTY E 360 - ISDA 2006 4.16(g) 30/360 ISMA, ICMA
* THIRTY E 360 ISDA - ISDA 2006 4.16(h)
* THIRTY E PLUS 360 - A month has 30 days. It rolls D2 to next month if D2 = 31
* ACT ACT ISDA - Splits accrued period into leap and non-leap year portions.
* ACT ACT ICMA - Used for US Treasury notes and bonds. Takes 3 dates and a annual_frequency.
* ACT 365 F - Denominator is always Fixed at 365, even in a leap year
* ACT 360 - Day difference divided by 360 - always
* ACT 365L - the 29 Feb is counted if it is in the date range

