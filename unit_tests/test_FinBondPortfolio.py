########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import datetime as dt

from financepy.utils.date import Date, from_datetime
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.bonds.bond import Bond

settlement = Date(19, 9, 2012)


def test_1():
    dc_type = DayCountTypes.THIRTY_360_BOND
    maturity_dt = Date(7, 3, 2013)
    coupon = 0.045
    clean_price = 101.99500000

    issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond.accrued_int, 4) == 0.001500
    assert round(ytm * 100, 4) == 0.2203


def test_2():
    dc_type = DayCountTypes.THIRTY_360_BOND
    maturity_dt = Date(7, 3, 2013)
    coupon = 0.045
    clean_price = 101.99500000

    issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond.accrued_int, 4) == 0.001500
    assert round(ytm * 100, 4) == 0.2203


def test_3():
    dc_type = DayCountTypes.THIRTY_E_360
    maturity_dt = Date(27, 9, 2013)
    coupon = 0.080000
    clean_price = 107.92000000

    issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond.accrued_int, 4) == 0.0382
    assert round(ytm * 100, 4) == 0.2380


def test_4():
    dc_type = DayCountTypes.THIRTY_E_360_ISDA
    maturity_dt = Date(7, 3, 2014)
    coupon = 0.022500
    clean_price = 102.9750

    issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond.accrued_int, 4) == 0.0008
    assert round(ytm * 100, 4) == 0.2172


def test_5():
    dc_type = DayCountTypes.THIRTY_E_PLUS_360
    maturity_dt = Date(7, 9, 2014)
    coupon = 0.0500000000
    clean_price = 109.35500000

    issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond.accrued_int, 4) == 0.0017
    assert round(ytm * 100, 4) == 0.2297


def test_6():
    dc_type = DayCountTypes.ACT_ACT_ISDA
    maturity_dt = Date(22, 1, 2015)
    coupon = 0.0275000000
    clean_price = 105.62500000

    issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond.accrued_int, 4) == 0.0044
    assert round(ytm * 100, 4) == 0.3334


def test_7():
    dc_type = DayCountTypes.ACT_ACT_ICMA
    maturity_dt = Date(7, 9, 2015)
    coupon = 0.0475000000
    clean_price = 112.98000000

    issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond.accrued_int, 4) == 0.0016
    assert round(ytm * 100, 4) == 0.3485


def test_8():
    dc_type = DayCountTypes.ACT_365F
    maturity_dt = Date(7, 12, 2015)
    coupon = 0.0800000000
    clean_price = 124.47000000

    issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond.accrued_int, 4) == 0.0228
    assert round(ytm * 100, 4) == 0.3405


def test_9():
    dc_type = DayCountTypes.ACT_360
    maturity_dt = Date(22, 1, 2016)
    coupon = 0.0200000000
    clean_price = 104.98000000

    issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond.accrued_int, 4) == 0.0033
    assert round(ytm * 100, 4) == 0.4930


def test_10():
    dc_type = DayCountTypes.ACT_365L
    maturity_dt = Date(7, 9, 2016)
    coupon = 0.0400000000
    clean_price = 113.49500000

    issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond.accrued_int, 4) == 0.0013
    assert round(ytm * 100, 4) == 0.5559


def test_11():
    dc_type = DayCountTypes.SIMPLE
    maturity_dt = Date(25, 8, 2017)
    coupon = 0.0875000000
    clean_price = 138.57000000

    issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond.accrued_int, 4) == 0.0060
    assert round(ytm * 100, 4) == 0.7652
