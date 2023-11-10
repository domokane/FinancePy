###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import os
import datetime as dt

from financepy.utils.date import Date, from_datetime
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.bonds.bond import Bond

settlement = Date(19, 9, 2012)


def test_1():
    accrual_type = DayCountTypes.THIRTY_360_BOND
    maturityDt = Date(7, 3, 2013)
    coupon = 0.045
    clean_price = 101.99500000

    issueDt = Date(maturityDt._d, maturityDt._m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issueDt, maturityDt,
                coupon, freq_type, accrual_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond._accrued_interest, 4) == 0.001500
    assert round(ytm * 100, 4) == 0.2203


def test_2():
    accrual_type = DayCountTypes.THIRTY_360_BOND
    maturityDt = Date(7, 3, 2013)
    coupon = 0.045
    clean_price = 101.99500000

    issueDt = Date(maturityDt._d, maturityDt._m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issueDt, maturityDt,
                coupon, freq_type, accrual_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond._accrued_interest, 4) == 0.001500
    assert round(ytm * 100, 4) == 0.2203


def test_3():
    accrual_type = DayCountTypes.THIRTY_E_360
    maturityDt = Date(27, 9, 2013)
    coupon = 0.080000
    clean_price = 107.92000000

    issueDt = Date(maturityDt._d, maturityDt._m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issueDt, maturityDt,
                coupon, freq_type, accrual_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond._accrued_interest, 4) == 0.0382
    assert round(ytm * 100, 4) == 0.2380


def test_4():
    accrual_type = DayCountTypes.THIRTY_E_360_ISDA
    maturityDt = Date(7, 3, 2014)
    coupon = 0.022500
    clean_price = 102.9750

    issueDt = Date(maturityDt._d, maturityDt._m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issueDt, maturityDt,
                coupon, freq_type, accrual_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond._accrued_interest, 4) == 0.0008
    assert round(ytm * 100, 4) == 0.2172


def test_5():
    accrual_type = DayCountTypes.THIRTY_E_PLUS_360
    maturityDt = Date(7, 9, 2014)
    coupon = 0.0500000000
    clean_price = 109.35500000

    issueDt = Date(maturityDt._d, maturityDt._m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issueDt, maturityDt,
                coupon, freq_type, accrual_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond._accrued_interest, 4) == 0.0017
    assert round(ytm * 100, 4) == 0.2297


def test_6():
    accrual_type = DayCountTypes.ACT_ACT_ISDA
    maturityDt = Date(22, 1, 2015)
    coupon = 0.0275000000
    clean_price = 105.62500000

    issueDt = Date(maturityDt._d, maturityDt._m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issueDt, maturityDt,
                coupon, freq_type, accrual_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond._accrued_interest, 4) == 0.0044
    assert round(ytm * 100, 4) == 0.3334


def test_7():
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    maturityDt = Date(7, 9, 2015)
    coupon = 0.0475000000
    clean_price = 112.98000000

    issueDt = Date(maturityDt._d, maturityDt._m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issueDt, maturityDt,
                coupon, freq_type, accrual_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond._accrued_interest, 4) == 0.0016
    assert round(ytm * 100, 4) == 0.3485


def test_8():
    accrual_type = DayCountTypes.ACT_365F
    maturityDt = Date(7, 12, 2015)
    coupon = 0.0800000000
    clean_price = 124.47000000

    issueDt = Date(maturityDt._d, maturityDt._m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issueDt, maturityDt,
                coupon, freq_type, accrual_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond._accrued_interest, 4) == 0.0228
    assert round(ytm * 100, 4) == 0.3405


def test_9():
    accrual_type = DayCountTypes.ACT_360
    maturityDt = Date(22, 1, 2016)
    coupon = 0.0200000000
    clean_price = 104.98000000

    issueDt = Date(maturityDt._d, maturityDt._m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issueDt, maturityDt,
                coupon, freq_type, accrual_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond._accrued_interest, 4) == 0.0033
    assert round(ytm * 100, 4) == 0.4930


def test_10():
    accrual_type = DayCountTypes.ACT_365L
    maturityDt = Date(7, 9, 2016)
    coupon = 0.0400000000
    clean_price = 113.49500000

    issueDt = Date(maturityDt._d, maturityDt._m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issueDt, maturityDt,
                coupon, freq_type, accrual_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond._accrued_interest, 4) == 0.0013
    assert round(ytm * 100, 4) == 0.5559


def test_11():
    accrual_type = DayCountTypes.SIMPLE
    maturityDt = Date(25, 8, 2017)
    coupon = 0.0875000000
    clean_price = 138.57000000

    issueDt = Date(maturityDt._d, maturityDt._m, 2000)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issueDt, maturityDt,
                coupon, freq_type, accrual_type)

    ytm = bond.yield_to_maturity(settlement, clean_price)
    assert round(bond._accrued_interest, 4) == 0.0060
    assert round(ytm * 100, 4) == 0.7652
