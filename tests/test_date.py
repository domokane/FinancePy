###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import datetime
import numpy as np
import time

from financepy.utils.date import Date, date_range

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
    assert Date(5, 1, 1900)._excel_date == 5
    assert Date(1, 3, 2020)._excel_date == 43891


# tests not refactored below
# - print() should be assert to value
# - do not need many values, just one call
# - some tests may be "parametrised" (probably not todo now)


def test_Date():

    startDate = Date(1, 1, 2018)

    # Checking for num_months = 0
    next_cds_date = startDate.next_cds_date(0)
    assert next_cds_date == Date(20, 3, 2018)

    # Checking for num_months = 120
    next_cds_date = startDate.next_cds_date(120)
    assert next_cds_date == Date(20, 3, 2028)

    startDate = Date(1, 1, 2018)

    # Checking for num_days = 1
    next_imm_date = startDate.add_days(1).next_imm_date()
    assert next_imm_date == Date(21, 3, 2018)

    # Checking for num_days = 365
    next_imm_date = startDate.add_days(365).next_imm_date()
    assert next_imm_date == Date(20, 3, 2019)


def test_DateTenors():

    startDate = Date(23, 2, 2018)

    tenor = "5d"
    assert startDate.add_tenor(tenor) == Date(28, 2, 2018)

    tenor = "7D"
    assert startDate.add_tenor(tenor) == Date(2, 3, 2018)

    tenor = "1W"
    assert startDate.add_tenor(tenor) == Date(2, 3, 2018)

    tenor = "4W"
    assert startDate.add_tenor(tenor) == Date(23, 3, 2018)

    tenor = "1M"
    assert startDate.add_tenor(tenor) == Date(23, 3, 2018)

    tenor = "24M"
    assert startDate.add_tenor(tenor) == Date(23, 2, 2020)

    tenor = "2Y"
    assert startDate.add_tenor(tenor) == Date(23, 2, 2020)

    tenor = "10y"
    assert startDate.add_tenor(tenor) == Date(23, 2, 2028)

    tenor = "0m"
    assert startDate.add_tenor(tenor) == Date(23, 2, 2018)

    tenor = "20Y"
    assert startDate.add_tenor(tenor) == Date(23, 2, 2038)

    tenor = "-5d"
    assert startDate.add_tenor(tenor) == Date(18, 2, 2018)

    tenor = "-7D"
    assert startDate.add_tenor(tenor) == Date(16, 2, 2018)

    tenor = "-1W"
    assert startDate.add_tenor(tenor) == Date(16, 2, 2018)

    tenor = "-4W"
    assert startDate.add_tenor(tenor) == Date(26, 1, 2018)

    tenor = "-1M"
    assert startDate.add_tenor(tenor) == Date(23, 1, 2018)

    tenor = "-24M"
    assert startDate.add_tenor(tenor) == Date(23, 2, 2016)

    tenor = "-2Y"
    assert startDate.add_tenor(tenor) == Date(23, 2, 2016)

    tenor = "-10y"
    assert startDate.add_tenor(tenor) == Date(23, 2, 2008)

    tenor = "-0m"
    assert startDate.add_tenor(tenor) == Date(23, 2, 2018)

    tenor = "-20Y"
    assert startDate.add_tenor(tenor) == Date(23, 2, 1998)


def test_DateRange():

    startDate = Date(1, 1, 2010)

    # Default
    endDate = startDate.add_days(3)
    dtRange = date_range(startDate, endDate)
    assert dtRange[0] == Date(1, 1, 2010)
    assert dtRange[-1] == Date(4, 1, 2010)

    # 1W Tenor
    endDate = startDate.add_days(20)
    tenor = "1W"
    dtRange = date_range(startDate, endDate, tenor)
    assert dtRange[0] == Date(1, 1, 2010)
    assert dtRange[-1] == Date(21, 1, 2010)

    # 7D Tenor
    tenor = "7D"
    dtRange = date_range(startDate, endDate, tenor)
    assert dtRange[1] == Date(8, 1, 2010)
    assert dtRange[2] == Date(15, 1, 2010)

    # Same startDate
    assert date_range(startDate, startDate)[0] == Date(1, 1, 2010)

    # startDate before endDate"
    assert len(date_range(endDate, startDate)) == 0


def test_DateAddMonths():

    startDate = Date(1, 1, 2010)

    months = [1, 3, 6, 9, 12, 24, 36, 48, 60]

    dates = startDate.add_months(months)

    assert dates[0] == Date(1, 2, 2010)
    assert dates[-1] == Date(1, 1, 2015)
    assert len(dates) == len(months)


def test_DateAddYears():

    startDate = Date(1, 1, 2010)

    # Simple list as input
    years = [1, 3, 5, 7, 10]
    dates_list = startDate.add_years(years)
    assert len(dates_list) == len(years)
    assert dates_list[0] == Date(1, 1, 2011)
    assert dates_list[-1] == Date(1, 1, 2020)

    # Numpy array as input
    years = np.array([1, 3, 5, 7, 10])
    dates_numpy = startDate.add_years(years)
    assert len(dates_numpy) == len(years)
    assert dates_numpy[0] == Date(1, 1, 2011)
    assert dates_numpy[-1] == Date(1, 1, 2020)

    # Fractional years as input
    years = np.array([1.5, 3.25, 5.75, 7.25, 10.0])
    dates_fractional = startDate.add_years(years)
    assert len(dates_fractional) == len(years)
    assert dates_fractional[0] == Date(1, 7, 2011)
    assert dates_fractional[-1] == Date(1, 1, 2020)

    # Fractional years with date foramtting in a numpy array as input
    dt = 1.0 / 365.0
    years = np.array(
        [1.5 + 2.0 * dt, 3.5 - 6 * dt, 5.75 + 3 * dt, 7.25 + dt, 10.0 + dt]
    )
    dates_fractional_np = startDate.add_years(years)
    assert len(dates_fractional_np) == len(years)
    assert dates_fractional_np[0] == Date(3, 7, 2011)
    assert dates_fractional_np[-1] == Date(2, 1, 2020)


def test_DateFormat():

    dt = Date(20, 10, 2019)

    # Date format in Bloomberg
    set_date_format(DateFormatTypes.BLOOMBERG)
    assert str(dt) == "10/20/19"

    # Date format in Datetime
    set_date_format(DateFormatTypes.DATETIME)
    assert str(dt) == "20/10/2019 00:00:00"

    # Date format in UK_LONGEST
    set_date_format(DateFormatTypes.UK_LONGEST)
    assert str(dt) == "SUN 20 OCT 2019"


def test_IntraDay():

    d1 = Date(20, 10, 2019, 0, 0, 0)
    d2 = Date(25, 10, 2019, 0, 0, 0)
    diff = d2 - d1
    assert round(diff, 4) == 5

    ###########################################################################

    d1 = Date(20, 10, 2019, 10, 0, 0)
    d2 = Date(25, 10, 2019, 10, 25, 0)
    diff = d2 - d1
    assert round(diff, 4) == 5.0174

    ###########################################################################

    d1 = Date(20, 10, 2019, 10, 0, 0)
    d2 = Date(20, 10, 2019, 10, 25, 30)
    diff = d2 - d1
    assert round(diff, 4) == 0.0177

    ###########################################################################

    d1 = Date(19, 10, 2019, 10, 0, 0)
    d2 = Date(20, 10, 2019, 10, 25, 40)
    diff = d2 - d1
    assert round(diff, 4) == 1.0178


def test_DateEOM():

    dt = Date(29, 2, 2000)
    assert dt.is_eom() is True

    dt = Date(28, 2, 2001)
    assert dt.is_eom() is True

    dt = Date(29, 2, 2004)
    assert dt.is_eom() is True

    dt = Date(28, 2, 2005)
    assert dt.is_eom() is True

    dt = Date(31, 3, 2003)
    assert dt.is_eom() is True

    dt = Date(30, 4, 2004)
    assert dt.is_eom() is True

    dt = Date(31, 5, 2004)
    assert dt.is_eom() is True

    dt = Date(31, 12, 2010)
    assert dt.is_eom() is True

    dt = Date(2, 2, 2000)
    assert dt.eom().is_eom() is True

    dt = Date(24, 2, 2001)
    assert dt.eom().is_eom() is True

    dt = Date(22, 2, 2004)
    assert dt.eom().is_eom() is True

    dt = Date(1, 2, 2005)
    assert dt.eom().is_eom() is True

    dt = Date(1, 3, 2003)
    assert dt.eom().is_eom() is True

    dt = Date(3, 4, 2004)
    assert dt.eom().is_eom() is True

    dt = Date(5, 5, 2004)
    assert dt.eom().is_eom() is True

    dt = Date(7, 12, 2010)
    assert dt.eom().is_eom() is True


def test_datetime():

    dt = Date(30, 12, 2021)
    assert dt.datetime()


def test_from_date():
    y, m, d = 2022, 11, 8
    dt1 = Date(d, m, y)
    dt2 = Date.from_date(datetime.date(y, m, d))
    assert dt1 == dt2


def test_list_of_date():
    dates = [Date(1, 1, 2020),
             Date(1, 2, 2020),
             Date(4, 3, 2020),
             Date(1, 1, 2021),
             Date(1, 1, 2022),
             ]

    # Test logical operations
    assert all(dates < Date(1, 1, 2023))
    assert all(dates <= Date(1, 1, 2022))
    assert all(dates > Date(1, 1, 1970))
    assert all(dates >= Date(1, 1, 2020))
    assert (dates >= Date(1, 1, 2021)) == [False, False, False, True, True]
    assert (dates == Date(1, 1, 2020)) == [True, False, False, False, False]

    # Test list of length 1
    assert ([Date(1, 1, 2020)] == Date(1, 1, 2020)) == [True]
    assert ([Date(1, 1, 2020)] == Date(1, 1, 1970)) == [False]

    # Test finding date difference
    assert (Date(1, 1, 2019) - dates) == [Date(1, 1, 2019) - d for d in dates]
    assert (dates - Date(1, 1, 2019)) == [Date(1, 1, 2019) - d for d in dates]
