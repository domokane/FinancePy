# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path
_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause
_install_double_click_pause()
import sys
import time
import datetime

import numpy as np

import add_fp_to_path

from financepy.utils.date_format import set_date_format
from financepy.utils.date_format import DateFormatTypes
from financepy.utils.date import Date, date_range
from financepy.utils import from_datetime




set_date_format(DateFormatTypes.UK_LONGEST)

########################################################################################


def test_dt():

    start_dt = Date(1, 1, 2018)

    assert Date(1, 1, 2018) == Date.from_string("1-1-2018", "%d-%m-%Y")

    print("DATE", "MONTHS", "CDS DATE")

    for num_months in range(0, 120):
        next_cds_date = start_dt.next_cds_date(num_months)
        print(str(start_dt), num_months, str(next_cds_date))

    start_dt = Date(1, 1, 2018)

    print("STARTDATE", "MONTHS", "CDS DATE")

    for num_months in range(0, 365):
        start_dt = start_dt.add_days(1)
        next_imm_date = start_dt.next_imm_date()
        print(num_months, str(start_dt), str(next_imm_date))


########################################################################################


def test_dt_tenors():

    start_dt = Date(23, 2, 2018)

    print("TENOR", "DATE")
    tenor = "5d"
    print(tenor, start_dt.add_tenor(tenor))

    tenor = "7D"
    print(tenor, start_dt.add_tenor(tenor))

    tenor = "1W"
    print(tenor, start_dt.add_tenor(tenor))

    tenor = "4W"
    print(tenor, start_dt.add_tenor(tenor))

    tenor = "1M"
    print(tenor, start_dt.add_tenor(tenor))

    tenor = "24M"
    print(tenor, start_dt.add_tenor(tenor))

    tenor = "2Y"
    print(tenor, start_dt.add_tenor(tenor))

    tenor = "10y"
    print(tenor, start_dt.add_tenor(tenor))

    tenor = "0m"
    print(tenor, start_dt.add_tenor(tenor))

    tenor = "20Y"
    print(tenor, start_dt.add_tenor(tenor))


########################################################################################


def test_dt_range():

    start_dt = Date(1, 1, 2010)

    print("Tenor", "Dates")

    end_dt = start_dt.add_days(3)
    tenor = "Default"
    print(tenor, date_range(start_dt, end_dt))

    end_dt = start_dt.add_days(20)

    tenor = "1W"
    print(tenor, date_range(start_dt, end_dt, tenor))

    tenor = "7D"
    print(tenor, date_range(start_dt, end_dt, tenor))

    print("Case", "Dates")

    case = "Same start_dt"
    print(case, date_range(start_dt, start_dt))
    case = "start_dt before end_dt"
    print(case, date_range(end_dt, start_dt))


########################################################################################


def test_dt_add_months():

    start_dt = Date(1, 1, 2010)

    print("Months", "Dates")

    months = [1, 3, 6, 9, 12, 24, 36, 48, 60]

    dates = start_dt.add_months(months)

    print("DATES", "DATE")

    for dt in dates:
        print("DATE", dt)


########################################################################################


def test_dt_add_years():

    start_dt = Date(1, 1, 2010)

    print("Years", "Dates")

    years = [1, 3, 5, 7, 10]
    dates1 = start_dt.add_years(years)
    for dt in dates1:
        print("DATES1", dt)

    years = np.array([1, 3, 5, 7, 10])
    dates2 = start_dt.add_years(years)
    for dt in dates2:
        print("DATES2", dt)

    years = np.array([1.5, 3.25, 5.75, 7.25, 10.0])
    dates3 = start_dt.add_years(years)

    for dt in dates3:
        print("DATES3", dt)

    dt = 1.0 / 365.0
    years = np.array(
        [1.5 + 2.0 * dt, 3.5 - 6 * dt, 5.75 + 3 * dt, 7.25 + dt, 10.0 + dt]
    )
    dates4 = start_dt.add_years(years)

    for dt in dates4:
        print("DATES4", dt)


########################################################################################


def test_dt_speed():

    num_steps = 100
    start = time.time()
    date_list = []
    for _ in range(0, num_steps):
        start_dt = Date(1, 1, 2010)
        date_list.append(start_dt)
    end = time.time()
    elapsed = end - start

    print("LABEL", "TIME")
    print("TIMING", elapsed)

    mem = sys.getsizeof(date_list)
    print("Mem:", mem)


########################################################################################


def test_dt_format():

    dt = Date(20, 10, 2019)
    print("FORMAT", "DATE")

    for format_type in DateFormatTypes:
        set_date_format(format_type)
        print(format_type.name, dt)


########################################################################################


def test_intra_day():

    print("Date1", "Date2", "Diff")
    d1 = Date(20, 10, 2019, 0, 0, 0)
    d2 = Date(25, 10, 2019, 0, 0, 0)
    diff = d2 - d1
    print(d1, d2, diff)
    print(d1.excel_dt, d2.excel_dt, diff)

    d1 = Date(20, 10, 2019, 10, 0, 0)
    d2 = Date(25, 10, 2019, 10, 25, 0)
    diff = d2 - d1
    print(d1, d2, diff)
    print(d1.excel_dt, d2.excel_dt, diff)

    d1 = Date(20, 10, 2019, 10, 0, 0)
    d2 = Date(20, 10, 2019, 10, 25, 30)
    diff = d2 - d1
    print(d1, d2, diff)
    print(d1.excel_dt, d2.excel_dt, diff)

    d1 = Date(19, 10, 2019, 10, 0, 0)
    d2 = Date(20, 10, 2019, 10, 25, 40)
    diff = d2 - d1
    print(d1, d2, diff)
    print(d1.excel_dt, d2.excel_dt, diff)


########################################################################################


def test_dt_eom():

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


########################################################################################


def test_add_weekdays():

    today = datetime.date(2022, 2, 13)  # Sunday 13th Feb
    next_weekday = from_datetime(today).add_weekdays(1)
    last_weekday = from_datetime(today).add_weekdays(-1)
    assert (last_weekday == Date(11, 2, 2022)) is True
    assert (next_weekday == Date(14, 2, 2022)) is True

    today = datetime.date(2022, 2, 13)  # Sunday 13th Feb
    next_weekday = from_datetime(today).add_weekdays(7)
    last_weekday = from_datetime(today).add_weekdays(-7)
    assert (last_weekday == Date(3, 2, 2022)) is True
    assert (next_weekday == Date(22, 2, 2022)) is True


########################################################################################

test_add_weekdays()

start = time.time()

test_dt()
test_dt_tenors()
test_dt_range()
test_dt_add_months()
test_dt_add_years()
test_dt_speed()
test_dt_format()
test_intra_day()
test_dt_eom()

end = time.time()
elapsed = end - start
# print("Elapsed time:", elapsed)


set_date_format(DateFormatTypes.UK_LONG)
