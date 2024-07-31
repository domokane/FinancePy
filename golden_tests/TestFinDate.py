###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
import time

import sys
sys.path.append("..")

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.date import set_date_format
from financepy.utils.date import DateFormatTypes
from financepy.utils.date import Date, date_range

test_cases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

set_date_format(DateFormatTypes.UK_LONGEST)


def test_dt():

    start_dt = Date(1, 1, 2018)

    assert Date(1, 1, 2018) == Date.from_string('1-1-2018', '%d-%m-%Y')

    test_cases.header("DATE", "MONTHS", "CDS DATE")

    for num_months in range(0, 120):
        next_cds_date = start_dt.next_cds_date(num_months)
        test_cases.print(str(start_dt), num_months, str(next_cds_date))

    start_dt = Date(1, 1, 2018)

    test_cases.header("STARTDATE", "MONTHS", "CDS DATE")

    for num_months in range(0, 365):
        start_dt = start_dt.add_days(1)
        next_imm_date = start_dt.next_imm_date()
        test_cases.print(num_months, str(start_dt), str(next_imm_date))

###############################################################################


def test_dtTenors():

    start_dt = Date(23, 2, 2018)

    test_cases.header("TENOR", "DATE")
    tenor = "5d"
    test_cases.print(tenor, start_dt.add_tenor(tenor))

    tenor = "7D"
    test_cases.print(tenor, start_dt.add_tenor(tenor))

    tenor = "1W"
    test_cases.print(tenor, start_dt.add_tenor(tenor))

    tenor = "4W"
    test_cases.print(tenor, start_dt.add_tenor(tenor))

    tenor = "1M"
    test_cases.print(tenor, start_dt.add_tenor(tenor))

    tenor = "24M"
    test_cases.print(tenor, start_dt.add_tenor(tenor))

    tenor = "2Y"
    test_cases.print(tenor, start_dt.add_tenor(tenor))

    tenor = "10y"
    test_cases.print(tenor, start_dt.add_tenor(tenor))

    tenor = "0m"
    test_cases.print(tenor, start_dt.add_tenor(tenor))

    tenor = "20Y"
    test_cases.print(tenor, start_dt.add_tenor(tenor))

###############################################################################


def test_dtRange():

    start_dt = Date(1, 1, 2010)

    test_cases.header("Tenor", "Dates")

    end_dt = start_dt.add_days(3)
    tenor = "Default"
    test_cases.print(tenor, date_range(start_dt, end_dt))

    end_dt = start_dt.add_days(20)
    tenor = "1W"
    test_cases.print(tenor, date_range(start_dt, end_dt, tenor))

    tenor = "7D"
    test_cases.print(tenor, date_range(start_dt, end_dt, tenor))

    test_cases.header("Case", "Dates")

    case = "Same start_dt"
    test_cases.print(case, date_range(start_dt, start_dt))
    case = "start_dt before end_dt"
    test_cases.print(case, date_range(end_dt, start_dt))

###############################################################################


def test_dtAddMonths():

    start_dt = Date(1, 1, 2010)

    test_cases.header("Months", "Dates")

    months = [1, 3, 6, 9, 12, 24, 36, 48, 60]

    dates = start_dt.add_months(months)

    test_cases.header("DATES", "DATE")

    for dt in dates:
        test_cases.print("DATE", dt)

###############################################################################


def test_dtAddYears():

    start_dt = Date(1, 1, 2010)

    test_cases.header("Years", "Dates")

    years = [1, 3, 5, 7, 10]
    dates1 = start_dt.add_years(years)
    for dt in dates1:
        test_cases.print("DATES1", dt)

    years = np.array([1, 3, 5, 7, 10])
    dates2 = start_dt.add_years(years)
    for dt in dates2:
        test_cases.print("DATES2", dt)

    years = np.array([1.5, 3.25, 5.75, 7.25, 10.0])
    dates3 = start_dt.add_years(years)

    for dt in dates3:
        test_cases.print("DATES3", dt)

    dt = 1.0/365.0
    years = np.array([1.5+2.0*dt, 3.5-6*dt, 5.75+3*dt, 7.25+dt, 10.0+dt])
    dates4 = start_dt.add_years(years)

    for dt in dates4:
        test_cases.print("DATES4", dt)

###############################################################################


def test_dtSpeed():

    num_steps = 100
    start = time.time()
    date_list = []
    for _ in range(0, num_steps):
        start_dt = Date(1, 1, 2010)
        date_list.append(start_dt)
    end = time.time()
    elapsed = end - start

    test_cases.header("LABEL", "TIME")
    test_cases.print("TIMING", elapsed)

    mem = sys.getsizeof(date_list)
    test_cases.print("Mem:", mem)


###############################################################################


def test_dtFormat():

    dt = Date(20, 10, 2019)
    test_cases.header("FORMAT", "DATE")

    for format_type in DateFormatTypes:
        set_date_format(format_type)
        test_cases.print(format_type.name, dt)

###############################################################################


def test_IntraDay():

    test_cases.header("Date1", "Date2", "Diff")
    d1 = Date(20, 10, 2019, 0, 0, 0)
    d2 = Date(25, 10, 2019, 0, 0, 0)
    diff = d2 - d1
    test_cases.print(d1, d2, diff)
    test_cases.print(d1.excel_dt, d2.excel_dt, diff)

    ###########################################################################

    d1 = Date(20, 10, 2019, 10, 0, 0)
    d2 = Date(25, 10, 2019, 10, 25, 0)
    diff = d2 - d1
    test_cases.print(d1, d2, diff)
    test_cases.print(d1.excel_dt, d2.excel_dt, diff)

    ###########################################################################

    d1 = Date(20, 10, 2019, 10, 0, 0)
    d2 = Date(20, 10, 2019, 10, 25, 30)
    diff = d2 - d1
    test_cases.print(d1, d2, diff)
    test_cases.print(d1.excel_dt, d2.excel_dt, diff)

    ###########################################################################

    d1 = Date(19, 10, 2019, 10, 0, 0)
    d2 = Date(20, 10, 2019, 10, 25, 40)
    diff = d2 - d1
    test_cases.print(d1, d2, diff)
    test_cases.print(d1.excel_dt, d2.excel_dt, diff)

###############################################################################


def test_dtEOM():

    dt = Date(29, 2, 2000)
    assert(dt.is_eom() is True)

    dt = Date(28, 2, 2001)
    assert(dt.is_eom() is True)

    dt = Date(29, 2, 2004)
    assert(dt.is_eom() is True)

    dt = Date(28, 2, 2005)
    assert(dt.is_eom() is True)

    dt = Date(31, 3, 2003)
    assert(dt.is_eom() is True)

    dt = Date(30, 4, 2004)
    assert(dt.is_eom() is True)

    dt = Date(31, 5, 2004)
    assert(dt.is_eom() is True)

    dt = Date(31, 12, 2010)
    assert(dt.is_eom() is True)

    dt = Date(2, 2, 2000)
    assert(dt.eom().is_eom() is True)

    dt = Date(24, 2, 2001)
    assert(dt.eom().is_eom() is True)

    dt = Date(22, 2, 2004)
    assert(dt.eom().is_eom() is True)

    dt = Date(1, 2, 2005)
    assert(dt.eom().is_eom() is True)

    dt = Date(1, 3, 2003)
    assert(dt.eom().is_eom() is True)

    dt = Date(3, 4, 2004)
    assert(dt.eom().is_eom() is True)

    dt = Date(5, 5, 2004)
    assert(dt.eom().is_eom() is True)

    dt = Date(7, 12, 2010)
    assert(dt.eom().is_eom() is True)

###############################################################################

import datetime
from financepy.utils import from_datetime

def test_add_weekdays():

    today = datetime.date(2022,2,13) # Sunday 13th Feb
    next_weekday = from_datetime(today).add_weekdays(1)
    last_weekday = from_datetime(today).add_weekdays(-1)
    assert( (last_weekday == Date(11, 2, 2022)) is True)
    assert( (next_weekday == Date(14, 2, 2022)) is True)

    today = datetime.date(2022,2,13) # Sunday 13th Feb
    next_weekday = from_datetime(today).add_weekdays(7)
    last_weekday = from_datetime(today).add_weekdays(-7)
    assert( (last_weekday == Date(3, 2, 2022)) is True)
    assert( (next_weekday == Date(22, 2, 2022)) is True)

###############################################################################

test_add_weekdays()

start = time.time()

test_dt()
test_dtTenors()
test_dtRange()
test_dtAddMonths()
test_dtAddYears()
test_dtSpeed()
test_dtFormat()
test_IntraDay()
test_dtEOM()

end = time.time()
elapsed = end - start
# print("Elapsed time:", elapsed)

test_cases.compareTestCases()

set_date_format(DateFormatTypes.UK_LONG)
