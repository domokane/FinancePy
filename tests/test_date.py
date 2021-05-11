###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
import time

from financepy.utils.date import Date, dateRange

# Not under test

from financepy.utils.date import DateFormatTypes, set_date_format

set_date_format(DateFormatTypes.UK_LONGEST)

# new sample tests


def test_add_days():
    assert Date(1, 1, 2018).add_days(-1).add_days(1) == Date(1, 1, 2018)


def test_from_string():
    assert Date.from_string("1-1-2018", "%d-%m-%Y") == Date(1, 1, 2018)


def test_weekday():
    assert Date(3, 3, 2021)._weekday == Date.WED


def test_excel_representation():
    assert Date(5, 1, 1900)._excelDate == 5
    assert Date(1, 3, 2020)._excelDate == 43891


# tests not refactored below
# - print() should be assert to value
# - do not need many values, just one call
# - some tests may be "parametrised" (probably not todo now)


def test_Date():

    startDate = Date(1, 1, 2018)

    for num_months in range(0, 120):
        next_cds_date = startDate.next_cds_date(num_months)
        print(str(startDate), num_months, str(next_cds_date))

    startDate = Date(1, 1, 2018)

    for num_months in range(0, 365):
        startDate = startDate.add_days(1)
        next_imm_date = startDate.next_imm_date()
        print(num_months, str(startDate), str(next_imm_date))


def test_DateTenors():

    startDate = Date(23, 2, 2018)

    tenor = "5d"
    print(tenor, startDate.add_tenor(tenor))

    tenor = "7D"
    print(tenor, startDate.add_tenor(tenor))

    tenor = "1W"
    print(tenor, startDate.add_tenor(tenor))

    tenor = "4W"
    print(tenor, startDate.add_tenor(tenor))

    tenor = "1M"
    print(tenor, startDate.add_tenor(tenor))

    tenor = "24M"
    print(tenor, startDate.add_tenor(tenor))

    tenor = "2Y"
    print(tenor, startDate.add_tenor(tenor))

    tenor = "10y"
    print(tenor, startDate.add_tenor(tenor))

    tenor = "0m"
    print(tenor, startDate.add_tenor(tenor))

    tenor = "20Y"
    print(tenor, startDate.add_tenor(tenor))


def test_DateRange():

    startDate = Date(1, 1, 2010)

    endDate = startDate.add_days(3)
    tenor = "Default"
    print(tenor, dateRange(startDate, endDate))

    endDate = startDate.add_days(20)
    tenor = "1W"
    print(tenor, dateRange(startDate, endDate, tenor))

    tenor = "7D"
    print(tenor, dateRange(startDate, endDate, tenor))

    case = "Same startDate"
    print(case, dateRange(startDate, startDate))
    case = "startDate before endDate"
    print(case, dateRange(endDate, startDate))


def test_DateAddMonths():

    startDate = Date(1, 1, 2010)

    months = [1, 3, 6, 9, 12, 24, 36, 48, 60]

    dates = startDate.add_months(months)

    for dt in dates:
        print("DATE", dt)


def test_DateAddYears():

    startDate = Date(1, 1, 2010)

    years = [1, 3, 5, 7, 10]
    dates1 = startDate.add_years(years)
    for dt in dates1:
        print("DATES1", dt)

    years = np.array([1, 3, 5, 7, 10])
    dates2 = startDate.add_years(years)
    for dt in dates2:
        print("DATES2", dt)

    years = np.array([1.5, 3.25, 5.75, 7.25, 10.0])
    dates3 = startDate.add_years(years)

    for dt in dates3:
        print("DATES3", dt)

    dt = 1.0 / 365.0
    years = np.array(
        [1.5 + 2.0 * dt, 3.5 - 6 * dt, 5.75 + 3 * dt, 7.25 + dt, 10.0 + dt]
    )
    dates4 = startDate.add_years(years)

    for dt in dates4:
        print("DATES4", dt)


def test_DateFormat():

    dt = Date(20, 10, 2019)

    for formatType in DateFormatTypes:
        set_date_format(formatType)
        print(formatType.name, dt)


def test_IntraDay():

    d1 = Date(20, 10, 2019, 0, 0, 0)
    d2 = Date(25, 10, 2019, 0, 0, 0)
    diff = d2 - d1
    print(d1, d2, diff)
    print(d1._excelDate, d2._excelDate, diff)

    ###########################################################################

    d1 = Date(20, 10, 2019, 10, 0, 0)
    d2 = Date(25, 10, 2019, 10, 25, 0)
    diff = d2 - d1
    print(d1, d2, diff)
    print(d1._excelDate, d2._excelDate, diff)

    ###########################################################################

    d1 = Date(20, 10, 2019, 10, 0, 0)
    d2 = Date(20, 10, 2019, 10, 25, 30)
    diff = d2 - d1
    print(d1, d2, diff)
    print(d1._excelDate, d2._excelDate, diff)

    ###########################################################################

    d1 = Date(19, 10, 2019, 10, 0, 0)
    d2 = Date(20, 10, 2019, 10, 25, 40)
    diff = d2 - d1
    print(d1, d2, diff)
    print(d1._excelDate, d2._excelDate, diff)


def test_DateEOM():

    dt = Date(29, 2, 2000)
    assert dt.isEOM() == True

    dt = Date(28, 2, 2001)
    assert dt.isEOM() == True

    dt = Date(29, 2, 2004)
    assert dt.isEOM() == True

    dt = Date(28, 2, 2005)
    assert dt.isEOM() == True

    dt = Date(31, 3, 2003)
    assert dt.isEOM() == True

    dt = Date(30, 4, 2004)
    assert dt.isEOM() == True

    dt = Date(31, 5, 2004)
    assert dt.isEOM() == True

    dt = Date(31, 12, 2010)
    assert dt.isEOM() == True

    dt = Date(2, 2, 2000)
    assert dt.EOM().isEOM() == True

    dt = Date(24, 2, 2001)
    assert dt.EOM().isEOM() == True

    dt = Date(22, 2, 2004)
    assert dt.EOM().isEOM() == True

    dt = Date(1, 2, 2005)
    assert dt.EOM().isEOM() == True

    dt = Date(1, 3, 2003)
    assert dt.EOM().isEOM() == True

    dt = Date(3, 4, 2004)
    assert dt.EOM().isEOM() == True

    dt = Date(5, 5, 2004)
    assert dt.EOM().isEOM() == True

    dt = Date(7, 12, 2010)
    assert dt.EOM().isEOM() == True

def test_datetime():
    
    dt = Date(30, 12, 2021)
    assert dt.datetime()
