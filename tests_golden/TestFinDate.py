###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.date import set_date_format
from financepy.utils.date import DateFormatTypes
from financepy.utils.date import Date, date_range
import numpy as np
import time

import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

set_date_format(DateFormatTypes.UK_LONGEST)


def test_Date():

    start_date = Date(1, 1, 2018)

    assert Date(1, 1, 2018) == Date.from_string('1-1-2018', '%d-%m-%Y')

    testCases.header("DATE", "MONTHS", "CDS DATE")

    for num_months in range(0, 120):
        next_cds_date = start_date.next_cds_date(num_months)
        testCases.print(str(start_date), num_months, str(next_cds_date))

    start_date = Date(1, 1, 2018)

    testCases.header("STARTDATE", "MONTHS", "CDS DATE")

    for num_months in range(0, 365):
        start_date = start_date.add_days(1)
        next_imm_date = start_date.next_imm_date()
        testCases.print(num_months, str(start_date), str(next_imm_date))

###############################################################################


def test_DateTenors():

    start_date = Date(23, 2, 2018)

    testCases.header("TENOR", "DATE")
    tenor = "5d"
    testCases.print(tenor, start_date.add_tenor(tenor))

    tenor = "7D"
    testCases.print(tenor, start_date.add_tenor(tenor))

    tenor = "1W"
    testCases.print(tenor, start_date.add_tenor(tenor))

    tenor = "4W"
    testCases.print(tenor, start_date.add_tenor(tenor))

    tenor = "1M"
    testCases.print(tenor, start_date.add_tenor(tenor))

    tenor = "24M"
    testCases.print(tenor, start_date.add_tenor(tenor))

    tenor = "2Y"
    testCases.print(tenor, start_date.add_tenor(tenor))

    tenor = "10y"
    testCases.print(tenor, start_date.add_tenor(tenor))

    tenor = "0m"
    testCases.print(tenor, start_date.add_tenor(tenor))

    tenor = "20Y"
    testCases.print(tenor, start_date.add_tenor(tenor))

###############################################################################


def test_DateRange():

    start_date = Date(1, 1, 2010)

    testCases.header("Tenor", "Dates")

    end_date = start_date.add_days(3)
    tenor = "Default"
    testCases.print(tenor, date_range(start_date, end_date))

    end_date = start_date.add_days(20)
    tenor = "1W"
    testCases.print(tenor, date_range(start_date, end_date, tenor))

    tenor = "7D"
    testCases.print(tenor, date_range(start_date, end_date, tenor))

    testCases.header("Case", "Dates")

    case = "Same start_date"
    testCases.print(case, date_range(start_date, start_date))
    case = "start_date before end_date"
    testCases.print(case, date_range(end_date, start_date))

###############################################################################


def test_DateAddMonths():

    start_date = Date(1, 1, 2010)

    testCases.header("Months", "Dates")

    months = [1, 3, 6, 9, 12, 24, 36, 48, 60]

    dates = start_date.add_months(months)

    testCases.header("DATES", "DATE")

    for dt in dates:
        testCases.print("DATE", dt)

###############################################################################


def test_DateAddYears():

    start_date = Date(1, 1, 2010)

    testCases.header("Years", "Dates")

    years = [1, 3, 5, 7, 10]
    dates1 = start_date.add_years(years)
    for dt in dates1:
        testCases.print("DATES1", dt)

    years = np.array([1, 3, 5, 7, 10])
    dates2 = start_date.add_years(years)
    for dt in dates2:
        testCases.print("DATES2", dt)

    years = np.array([1.5, 3.25, 5.75, 7.25, 10.0])
    dates3 = start_date.add_years(years)

    for dt in dates3:
        testCases.print("DATES3", dt)

    dt = 1.0/365.0
    years = np.array([1.5+2.0*dt, 3.5-6*dt, 5.75+3*dt, 7.25+dt, 10.0+dt])
    dates4 = start_date.add_years(years)

    for dt in dates4:
        testCases.print("DATES4", dt)

###############################################################################


def test_DateSpeed():

    num_steps = 100
    start = time.time()
    dateList = []
    for _ in range(0, num_steps):
        start_date = Date(1, 1, 2010)
        dateList.append(start_date)
    end = time.time()
    elapsed = end - start

    testCases.header("LABEL", "TIME")
    testCases.print("TIMING", elapsed)

    mem = sys.getsizeof(dateList)
    testCases.print("Mem:", mem)


###############################################################################


def test_DateFormat():

    dt = Date(20, 10, 2019)
    testCases.header("FORMAT", "DATE")

    for format_type in DateFormatTypes:
        set_date_format(format_type)
        testCases.print(format_type.name, dt)

###############################################################################


def test_IntraDay():

    testCases.header("Date1", "Date2", "Diff")
    d1 = Date(20, 10, 2019, 0, 0, 0)
    d2 = Date(25, 10, 2019, 0, 0, 0)
    diff = d2 - d1
    testCases.print(d1, d2, diff)
    testCases.print(d1._excel_date, d2._excel_date, diff)

    ###########################################################################

    d1 = Date(20, 10, 2019, 10, 0, 0)
    d2 = Date(25, 10, 2019, 10, 25, 0)
    diff = d2 - d1
    testCases.print(d1, d2, diff)
    testCases.print(d1._excel_date, d2._excel_date, diff)

    ###########################################################################

    d1 = Date(20, 10, 2019, 10, 0, 0)
    d2 = Date(20, 10, 2019, 10, 25, 30)
    diff = d2 - d1
    testCases.print(d1, d2, diff)
    testCases.print(d1._excel_date, d2._excel_date, diff)

    ###########################################################################

    d1 = Date(19, 10, 2019, 10, 0, 0)
    d2 = Date(20, 10, 2019, 10, 25, 40)
    diff = d2 - d1
    testCases.print(d1, d2, diff)
    testCases.print(d1._excel_date, d2._excel_date, diff)

###############################################################################


def test_DateEOM():

    dt = Date(29, 2, 2000)
    assert(dt.is_eom() == True)

    dt = Date(28, 2, 2001)
    assert(dt.is_eom() == True)

    dt = Date(29, 2, 2004)
    assert(dt.is_eom() == True)

    dt = Date(28, 2, 2005)
    assert(dt.is_eom() == True)

    dt = Date(31, 3, 2003)
    assert(dt.is_eom() == True)

    dt = Date(30, 4, 2004)
    assert(dt.is_eom() == True)

    dt = Date(31, 5, 2004)
    assert(dt.is_eom() == True)

    dt = Date(31, 12, 2010)
    assert(dt.is_eom() == True)

    dt = Date(2, 2, 2000)
    assert(dt.eom().is_eom() == True)

    dt = Date(24, 2, 2001)
    assert(dt.eom().is_eom() == True)

    dt = Date(22, 2, 2004)
    assert(dt.eom().is_eom() == True)

    dt = Date(1, 2, 2005)
    assert(dt.eom().is_eom() == True)

    dt = Date(1, 3, 2003)
    assert(dt.eom().is_eom() == True)

    dt = Date(3, 4, 2004)
    assert(dt.eom().is_eom() == True)

    dt = Date(5, 5, 2004)
    assert(dt.eom().is_eom() == True)

    dt = Date(7, 12, 2010)
    assert(dt.eom().is_eom() == True)

###############################################################################


start = time.time()

test_Date()
test_DateTenors()
test_DateRange()
test_DateAddMonths()
test_DateAddYears()
test_DateSpeed()
test_DateFormat()
test_IntraDay()
test_DateEOM()

end = time.time()
elapsed = end - start
# print("Elapsed time:", elapsed)

testCases.compareTestCases()

set_date_format(DateFormatTypes.UK_LONG)
