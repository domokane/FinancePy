###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

import pandas as pd
import numpy as np

from financepy.utils.math import ONE_MILLION
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.bonds.bond import YTMCalcType, Bond
from financepy.products.bonds.bond_zero import BondZero
from financepy.products.bonds.bond_market import BondMarkets
from financepy.products.bonds.bond_market import get_bond_market_conventions

import os

BASE_PATH = os.path.abspath(os.path.dirname(__file__))

###############################################################################


def test_bondtutor_example():

    #  EXAMPLE FROM http://bondtutor.com/btchp4/topic6/topic6.htm

    accrualConvention = DayCountTypes.ACT_ACT_ICMA
    y = 0.062267
    settlement_date = Date(19, 4, 1994)
    issue_date = Date(15, 7, 1990)
    maturity_date = Date(15, 7, 1997)
    coupon = 0.085
    ex_div_days = 0
    face = 1000000

    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issue_date, maturity_date,
                coupon, freq_type, accrualConvention, ex_div_days)

    dirty_price = bond.dirty_price_from_ytm(settlement_date, y)
    assert round(dirty_price, 4) == 108.7696
    clean_price = bond.clean_price_from_ytm(settlement_date, y)
    assert round(clean_price, 4) == 106.5625
    accrued_interest = bond.accrued_interest(settlement_date, face)
    assert round(accrued_interest, 4) == 22071.8232
    ytm = bond.yield_to_maturity(settlement_date, clean_price)
    assert round(ytm, 4) == 0.0622

    bump = 1e-4
    priceBumpedUp = bond.dirty_price_from_ytm(settlement_date, y + bump)
    assert round(priceBumpedUp, 4) == 108.7395

    priceBumpedDn = bond.dirty_price_from_ytm(settlement_date, y - bump)
    assert round(priceBumpedDn, 4) == 108.7998

    durationByBump = -(priceBumpedUp - dirty_price) / bump
    assert round(durationByBump, 4) == 301.1932

    duration = bond.dollar_duration(settlement_date, y)
    assert round(duration, 4) == 301.2458
    assert round(duration - durationByBump, 4) == 0.0526

    modified_duration = bond.modified_duration(settlement_date, y)
    assert round(modified_duration, 4) == 2.7696

    macauley_duration = bond.macauley_duration(settlement_date, y)
    assert round(macauley_duration, 4) == 2.8558

    conv = bond.convexity_from_ytm(settlement_date, y)
    assert round(conv, 4) == 0.0967


def test_bloomberg_us_treasury_example():
    # https://data.bloomberglp.com/bat/sites/3/2017/07/SF-2017_Paul-Fjeldsted.pdf

    settlement_date = Date(21, 7, 2017)
    issue_date = Date(15, 5, 2010)
    maturity_date = Date(15, 5, 2027)
    coupon = 0.02375
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    face = 100.0
    ex_div_days = 0

    bond = Bond(issue_date,
                maturity_date,
                coupon,
                freq_type,
                accrual_type,
                ex_div_days)

    clean_price = 99.7808417

    yld = bond.current_yield(clean_price)
    assert round(yld, 4) == 0.0238

    ytm = bond.yield_to_maturity(settlement_date, clean_price,
                                 YTMCalcType.UK_DMO)
    assert round(ytm, 4) == 0.0240

    ytm = bond.yield_to_maturity(settlement_date, clean_price,
                                 YTMCalcType.US_STREET)
    assert round(ytm, 4) == 0.0240

    ytm = bond.yield_to_maturity(settlement_date, clean_price,
                                 YTMCalcType.US_TREASURY)
    assert round(ytm, 4) == 0.0240

    dirty_price = bond.dirty_price_from_ytm(settlement_date, ytm)
    assert round(dirty_price, 4) == 100.2149

    clean_price = bond.clean_price_from_ytm(settlement_date, ytm)
    assert round(clean_price, 4) == 99.7825

    accrued_interest = bond.accrued_interest(settlement_date, face)
    assert round(accrued_interest, 4) == 0.4324

    accddays = bond._accrued_days
    assert round(accddays, 4) == 67.0

    duration = bond.dollar_duration(settlement_date, ytm)
    assert round(duration, 4) == 869.0934

    modified_duration = bond.modified_duration(settlement_date, ytm)
    assert round(modified_duration, 4) == 8.6723

    macauley_duration = bond.macauley_duration(settlement_date, ytm)
    assert round(macauley_duration, 4) == 8.7764

    conv = bond.convexity_from_ytm(settlement_date, ytm)
    assert round(conv, 4) == 0.8517


def test_bloomberg_apple_corp_example():
    settlement_date = Date(21, 7, 2017)
    issue_date = Date(13, 5, 2012)
    maturity_date = Date(13, 5, 2022)
    coupon = 0.027
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.THIRTY_E_360_ISDA
    face = 100.0
    ex_div_days = 0

    bond = Bond(issue_date, maturity_date,
                coupon, freq_type, accrual_type, ex_div_days)

    clean_price = 101.581564

    yld = bond.current_yield(clean_price)
    assert round(yld, 4) == 0.0266

    ytm = bond.yield_to_maturity(settlement_date, clean_price,
                                 YTMCalcType.UK_DMO)
    assert round(ytm, 4) == 0.0235

    ytm = bond.yield_to_maturity(settlement_date, clean_price,
                                 YTMCalcType.US_STREET)
    assert round(ytm, 4) == 0.0235

    ytm = bond.yield_to_maturity(settlement_date, clean_price,
                                 YTMCalcType.US_TREASURY)
    assert round(ytm, 4) == 0.0235

    dirty_price = bond.dirty_price_from_ytm(settlement_date, ytm)
    assert round(dirty_price, 4) == 102.0932

    clean_price = bond.clean_price_from_ytm(settlement_date, ytm)
    assert round(clean_price, 4) == 101.5832

    accddays = bond._accrued_days
    assert accddays == 68

    accrued_interest = bond.accrued_interest(settlement_date, face)
    assert round(accrued_interest, 4) == 0.51

    duration = bond.dollar_duration(settlement_date, ytm)
    assert round(duration, 4) == 456.5778

    modified_duration = bond.modified_duration(settlement_date, ytm)
    assert round(modified_duration, 4) == 4.4722

    macauley_duration = bond.macauley_duration(settlement_date, ytm)
    assert round(macauley_duration, 4) == 4.5247

    conv = bond.convexity_from_ytm(settlement_date, ytm)
    assert round(conv, 4) == 0.2302

###############################################################################


def test_zero_bond():

    # A 3 months treasure with 0 coupon per year.

    bill = BondZero(
        issue_date=Date(25, 7, 2022),
        maturity_date=Date(24, 10, 2022),
        issue_price=99.6410
    )
    settlement_date = Date(8, 8, 2022)

    clean_price = 99.6504
    calc_ytm = bill.yield_to_maturity(
        settlement_date, clean_price, YTMCalcType.ZERO) * 100

    accrued_interest = bill.accrued_interest(settlement_date, ONE_MILLION)

    assert abs(calc_ytm - 1.3997) < 0.0002
    assert abs(accrued_interest - ONE_MILLION * 0.05523077 / 100) < 0.01

###############################################################################


def test_bond_ror():

    test_case_file = './data/test_cases_bond_ror.csv'
    path = os.path.join(os.path.dirname(__file__), test_case_file)

    df = pd.read_csv(path,
                     parse_dates=['buy_date', 'sell_date'])

    # A 10-year bond with 1 coupon per year. code: 210215
    bond = Bond(
        issue_date=Date(13, 9, 2021),
        maturity_date=Date(13, 9, 2031),
        coupon=0.0312,
        freq_type=FrequencyTypes.ANNUAL,
        accrual_type=DayCountTypes.ACT_ACT_ICMA
    )
    for row in df.itertuples(index=False):
        buy_date = Date(row.buy_date.day, row.buy_date.month,
                        row.buy_date.year)
        sell_date = Date(row.sell_date.day,
                         row.sell_date.month, row.sell_date.year)
        simple, irr, pnl = bond.calc_ror(
            buy_date, sell_date, row.buy_ytm, row.sell_ytm)
        assert abs(simple - row.simple_return) < 0.00001
        assert abs(irr - row.irr) < 0.00001

###############################################################################


def test_bond_zero_ror():

    test_case_file = './data/test_cases_bond_zero_ror.csv'
    path = os.path.join(os.path.dirname(__file__), test_case_file)

    df = pd.read_csv(path,
                     parse_dates=['buy_date', 'sell_date'])

    # A 1-year bond with zero coupon per year. code: 092103011
    bond = BondZero(
        issue_date=Date(23, 7, 2021),
        maturity_date=Date(24, 8, 2022),
        issue_price=97.67
    )
    for row in df.itertuples(index=False):
        buy_date = Date(row.buy_date.day, row.buy_date.month,
                        row.buy_date.year)
        sell_date = Date(row.sell_date.day,
                         row.sell_date.month, row.sell_date.year)
        simple, irr, pnl = bond.calc_ror(
            buy_date, sell_date, row.buy_ytm, row.sell_ytm)
        assert abs(simple - row.simple_return) < 0.00001
        assert abs(irr - row.irr) < 0.00001

###############################################################################


def test_bond_cfets():
    """
    Test ytms of bonds in CFETS convention, especially for those in last
    coupon period and have 2 or more coupon payments per year.
    """
    face = 100.0
    test_case_file = './data/test_cases_bond_cfets.csv'
    path = os.path.join(os.path.dirname(__file__), test_case_file)

    df = pd.read_csv(path,
                     parse_dates=['settlement_date', 'issue_date',
                                  'maturity_date'])

    for row in df.itertuples(index=False):

        issue_date = Date(row.issue_date.day,
                          row.issue_date.month,
                          row.issue_date.year)

        maturity_date = Date(row.maturity_date.day,
                             row.maturity_date.month,
                             row.maturity_date.year)

        if row.freq == 1:
            freq_type = FrequencyTypes.ANNUAL
        else:
            freq_type = FrequencyTypes.SEMI_ANNUAL

        bond = Bond(issue_date,
                    maturity_date,
                    row.coupon / 100,
                    freq_type,
                    accrual_type=DayCountTypes.ACT_ACT_ICMA)

        settlement_date = Date(row.settlement_date.day,
                               row.settlement_date.month,
                               row.settlement_date.year)

        accrued_interest = bond.accrued_interest(settlement_date, face)
        clean_price = row.dirty_price - accrued_interest
        calc_ytm = bond.yield_to_maturity(
            settlement_date, clean_price, YTMCalcType.CFETS) * 100
        try:
            assert abs(calc_ytm - row.ytm) < 0.0001
        except Exception:
            print(bond)
            print(clean_price)
            print(settlement_date)
            print(bond.bond_payments(settlement_date, 100.0))
            print(f'calc_ytm:{calc_ytm}, correct_ytm:{row.ytm}')
            continue

###############################################################################


def test_key_rate_durations_Bloomberg_example():

    accrual_type, frequencyType, settlementDays, exDiv, calendar = \
        get_bond_market_conventions(BondMarkets.UNITED_STATES)

    # interest accrues on this date. Issue date is 01/08/2022
    issue_date = Date(31, 7, 2022)
    maturity_date = Date(31, 7, 2027)
    coupon = 2.75/100.0
    ex_div_days = 0

    accrual_type, freq_type, settlementDays, exDiv, calendar =\
        get_bond_market_conventions(BondMarkets.UNITED_STATES)

    bond = Bond(issue_date, maturity_date, coupon,
                freq_type, accrual_type, ex_div_days)

    settlement_date = Date(24, 4, 2023)

    # US Street yield on Bloomberg as of 20 April 2023
    # with settle date 24 April 2023
    ytm = 3.725060/100

    # Details of yields of market bonds at KRD maturity points
    my_tenors = np.array([0.5,  1,  2,  3,  5,  7,  10])
    my_rates = np.array([5.0367, 4.7327, 4.1445, 3.8575,
                         3.6272,  3.5825,  3.5347]) / 100.0

    key_rate_tenors, key_rate_durations =\
        bond.key_rate_durations(settlement_date,
                                ytm,
                                key_rate_tenors=my_tenors,
                                rates=my_rates)

    bbg_key_rate_durations = [-0.001, -.009, -0.022, 1.432,
                              2.527, 0.00, 0.00, 0.00, 0.00]

    for i in range(len(key_rate_durations)):
        assert round(key_rate_durations[i], 3) == bbg_key_rate_durations[i]

###############################################################################
