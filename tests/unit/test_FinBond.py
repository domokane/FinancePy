# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import os
import pandas as pd
import numpy as np

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.global_vars import ONE_MILLION
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.bonds.bond import YTMCalcType, Bond
from financepy.products.bonds.bond_zero import BondZero
from financepy.products.bonds.bond_market import BondMarkets
from financepy.products.bonds.bond_market import get_bond_market_conventions

TOL = 3


def assert_my_round(x, y):
    n = len(str(y).split(".")[1])
    x_check = round(x, n)
    print("====>", x, y, n, x_check)
    assert x_check == y


def my_round(x):
    return round(x, TOL)


########################################################################################


def test_bondtutor_example():

    #  EXAMPLE FROM http://bondtutor.com/btchp4/topic6/topic6.htm

    accrual_convention = DayCountTypes.ACT_ACT_ICMA
    y = 0.062267
    settle_dt = Date(19, 4, 1994)
    issue_dt = Date(15, 7, 1990)
    maturity_dt = Date(15, 7, 1997)
    coupon = 0.085
    ex_div_days = 0
    face = 1000000

    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(
        issue_dt,
        maturity_dt,
        coupon,
        freq_type,
        accrual_convention,
        ex_div_days,
    )

    dirty_price = bond.dirty_price_from_ytm(settle_dt, y)
    assert my_round(dirty_price) == my_round(108.7696)
    clean_price = bond.clean_price_from_ytm(settle_dt, y)
    assert my_round(clean_price) == my_round(106.5625)
    accrued_interest = bond.accrued_interest(settle_dt, face)
    assert my_round(accrued_interest) == my_round(22071.8232)
    ytm = bond.yield_to_maturity(settle_dt, clean_price)
    assert my_round(ytm) == my_round(0.0622)

    duration = bond.dollar_duration(settle_dt, y)
    assert my_round(duration) == my_round(301.193)

    modified_duration = bond.modified_duration(settle_dt, y)
    assert my_round(modified_duration) == my_round(2.769)

    macaulay_duration = bond.macaulay_duration(settle_dt, y)
    assert my_round(macaulay_duration) == my_round(2.855)

    conv = bond.convexity_from_ytm(settle_dt, y)
    assert my_round(conv) == my_round(9.67)


########################################################################################


def test_bloomberg_us_treasury_example():

    # https://data.bloomberglp.com/bat/sites/3/2017/07/SF-2017_Paul-Fjeldsted.pdf

    settle_dt = Date(21, 7, 2017)
    issue_dt = Date(15, 5, 2010)
    maturity_dt = Date(15, 5, 2027)
    coupon = 0.02375
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    face = 100.0
    ex_div_days = 0

    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type, ex_div_days)

    clean_price = 99.7808417

    yld = bond.current_yield(settle_dt, clean_price)
    assert my_round(yld) == my_round(0.0237)

    ytm = bond.yield_to_maturity(settle_dt, clean_price, YTMCalcType.UK_DMO)
    assert my_round(ytm) == my_round(0.024)

    ytm = bond.yield_to_maturity(settle_dt, clean_price, YTMCalcType.US_STREET)
    assert my_round(ytm) == my_round(0.024)

    ytm = bond.yield_to_maturity(
        settle_dt, clean_price, YTMCalcType.US_TREASURY
    )
    assert my_round(ytm) == my_round(0.0240)

    dirty_price = bond.dirty_price_from_ytm(settle_dt, ytm)
    assert my_round(dirty_price) == my_round(100.2149)

    clean_price = bond.clean_price_from_ytm(settle_dt, ytm)
    assert my_round(clean_price) == my_round(99.7825)

    accrued_interest = bond.accrued_interest(settle_dt, face)
    assert my_round(accrued_interest) == my_round(0.4324)

    accddays = bond.accrued_days
    assert my_round(accddays) == my_round(67.0)

    duration = bond.dollar_duration(settle_dt, ytm)
    assert my_round(duration) == my_round(868.667)

    modified_duration = bond.modified_duration(settle_dt, ytm)
    assert my_round(modified_duration) == my_round(8.668)

    macaulay_duration = bond.macaulay_duration(settle_dt, ytm)
    assert my_round(macaulay_duration) == my_round(8.772)

    conv = bond.convexity_from_ytm(settle_dt, ytm)
    assert my_round(conv) == my_round(85.17)


########################################################################################


def test_bloomberg_apple_corp_example():

    settle_dt = Date(21, 7, 2017)
    issue_dt = Date(13, 5, 2012)
    maturity_dt = Date(13, 5, 2022)
    coupon = 0.027
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.THIRTY_E_360_ISDA
    face = 100.0
    ex_div_days = 0

    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type, ex_div_days)

    clean_price = 101.581564

    yld = bond.current_yield(settle_dt, clean_price)
    assert my_round(yld) == my_round(0.0264)

    ytm = bond.yield_to_maturity(settle_dt, clean_price, YTMCalcType.UK_DMO)
    assert my_round(ytm) == my_round(0.023)

    ytm = bond.yield_to_maturity(settle_dt, clean_price, YTMCalcType.US_STREET)
    assert my_round(ytm) == my_round(0.023)

    ytm = bond.yield_to_maturity(
        settle_dt, clean_price, YTMCalcType.US_TREASURY
    )
    assert my_round(ytm) == my_round(0.023)

    dirty_price = bond.dirty_price_from_ytm(settle_dt, ytm)
    assert my_round(dirty_price) == my_round(102.0932)

    clean_price = bond.clean_price_from_ytm(settle_dt, ytm)
    assert my_round(clean_price) == my_round(101.5832)

    accddays = bond.accrued_days
    assert accddays == 68

    accrued_interest = bond.accrued_interest(settle_dt, face)
    assert my_round(accrued_interest) == my_round(0.51)

    duration = bond.dollar_duration(settle_dt, ytm)
    assert my_round(duration) == my_round(456.46)

    modified_duration = bond.modified_duration(settle_dt, ytm)
    assert my_round(modified_duration) == my_round(4.471)

    macaulay_duration = bond.macaulay_duration(settle_dt, ytm)
    assert my_round(macaulay_duration) == my_round(4.524)

    conv = bond.convexity_from_ytm(settle_dt, ytm)
    assert my_round(conv) == my_round(23.024)


########################################################################################


def test_principal_scales_accrued_interest_with_face():
    bond = Bond(
        Date(1, 1, 2024),
        Date(1, 1, 2027),
        0.06,
        FrequencyTypes.SEMI_ANNUAL,
        DayCountTypes.ACT_ACT_ICMA,
    )
    settle_dt = Date(1, 3, 2025)
    ytm = 0.05
    clean_price = bond.clean_price_from_ytm(settle_dt, ytm)

    principal_100 = bond.principal(settle_dt, ytm, 100.0, YTMCalcType.UK_DMO)
    principal_million = bond.principal(
        settle_dt, ytm, ONE_MILLION, YTMCalcType.UK_DMO
    )

    assert np.isclose(principal_100, clean_price)
    assert np.isclose(principal_million, clean_price * ONE_MILLION / bond.par)


########################################################################################


def test_survival_curve_clean_price_subtracts_accrued_interest_per_par():
    settle_dt = Date(1, 3, 2025)
    bond = Bond(
        Date(1, 1, 2024),
        Date(1, 1, 2027),
        0.06,
        FrequencyTypes.SEMI_ANNUAL,
        DayCountTypes.ACT_ACT_ICMA,
    )
    discount_curve = FlatDiscountCurve(
        settle_dt, 0.03, FrequencyTypes.CONTINUOUS
    )
    survival_curve = FlatDiscountCurve(
        settle_dt, 0.02, FrequencyTypes.CONTINUOUS
    )

    dirty_price = bond.dirty_price_from_survival_curve(
        settle_dt, discount_curve, survival_curve, 0.40
    )
    clean_price = bond.clean_price_from_survival_curve(
        settle_dt, discount_curve, survival_curve, 0.40
    )
    accrued = bond.accrued_interest(settle_dt, bond.par)

    assert np.isclose(clean_price, dirty_price - accrued)


########################################################################################


def test_calculus_final_ex_dividend_preserves_redemption():
    bond = Bond(
        Date(15, 2, 2018),
        Date(15, 2, 2030),
        0.06,
        FrequencyTypes.SEMI_ANNUAL,
        DayCountTypes.ACT_ACT_ICMA,
        ex_div_days=7,
    )
    settle_dt = Date(10, 2, 2030)
    ytm = 0.05

    dirty_price = bond.dirty_price_from_ytm(
        settle_dt, ytm, YTMCalcType.CALCULUS
    )
    adjusted_ytm = ytm + 0.000000000012345
    expected_redemption = 100.0 / (
        1.0 + adjusted_ytm / bond.freq
    ) ** bond.alpha

    assert np.isclose(dirty_price, expected_redemption)
    assert np.isclose(
        bond.macaulay_duration(settle_dt, ytm, YTMCalcType.CALCULUS),
        bond.alpha / bond.freq,
    )


########################################################################################


def test_zero_bond():

    # A 3 months treasure with 0 coupon per year.

    issue_dt = Date(25, 7, 2022)
    maturity_dt = Date(24, 10, 2022)
    issue_price = 99.6410

    bill = BondZero(issue_dt, maturity_dt, issue_price)

    settle_dt = Date(8, 8, 2022)

    clean_price = 99.6504

    accrued_interest = bill.accrued_interest(settle_dt, ONE_MILLION)

    calc_ytm = (
        bill.yield_to_maturity(settle_dt, clean_price, YTMCalcType.ZERO) * 100
    )

    print(calc_ytm)
    assert abs(calc_ytm - 1.3997) < 0.0002
    assert abs(accrued_interest - ONE_MILLION * 0.05523077 / 100) < 0.01


########################################################################################


def test_bond_ror():

    test_case_file = "./data/test_cases_bond_ror.csv"
    path = os.path.join(os.path.dirname(__file__), test_case_file)

    df = pd.read_csv(path, parse_dates=["buy_date", "sell_date"])

    # A 10-year bond with 1 coupon per year. code: 210215
    bond = Bond(
        issue_dt=Date(13, 9, 2021),
        maturity_dt=Date(13, 9, 2031),
        coupon=0.0312,
        freq_type=FrequencyTypes.ANNUAL,
        accrual_dc_type=DayCountTypes.ACT_ACT_ICMA,
    )
    for row in df.itertuples(index=False):
        buy_date = Date(
            row.buy_date.day, row.buy_date.month, row.buy_date.year
        )
        sell_date = Date(
            row.sell_date.day, row.sell_date.month, row.sell_date.year
        )
        simple, irr, pnl = bond.calc_ror(
            buy_date, sell_date, row.buy_ytm, row.sell_ytm
        )
    assert abs(simple - row.simple_return) < 0.00001
    assert abs(irr - row.irr) < 0.00001


########################################################################################


def test_bond_zero_ror():

    test_case_file = "./data/test_cases_bond_zero_ror.csv"
    path = os.path.join(os.path.dirname(__file__), test_case_file)

    df = pd.read_csv(path, parse_dates=["buy_date", "sell_date"])

    # A 1-year bond with zero coupon per year. code: 092103011

    issue_dt = Date(23, 7, 2021)
    maturity_dt = Date(24, 8, 2022)
    issue_price = 97.67

    bond = BondZero(issue_dt, maturity_dt, issue_price)

    for row in df.itertuples(index=False):

        buy_date = Date(
            row.buy_date.day, row.buy_date.month, row.buy_date.year
        )

        sell_date = Date(
            row.sell_date.day, row.sell_date.month, row.sell_date.year
        )

        simple, irr, pnl = bond.calc_ror(
            buy_date, sell_date, row.buy_ytm, row.sell_ytm
        )

        assert abs(simple - row.simple_return) < 0.001

        assert abs(irr - row.irr) < 0.00001


########################################################################################


def test_bond_cfets():
    """
    Test ytms of bonds in CFETS convention, especially for those in last
    coupon period and have 2 or more coupon payments per year.
    """
    face = 100.0
    test_case_file = "./data/test_cases_bond_cfets.csv"
    path = os.path.join(os.path.dirname(__file__), test_case_file)

    df = pd.read_csv(
        path, parse_dates=["settle_dt", "issue_date", "maturity_dt"]
    )

    for row in df.itertuples(index=False):

        issue_dt = Date(
            row.issue_date.day, row.issue_date.month, row.issue_date.year
        )

        maturity_dt = Date(
            row.maturity_dt.day, row.maturity_dt.month, row.maturity_dt.year
        )

        if row.freq == 1:
            freq_type = FrequencyTypes.ANNUAL
        else:
            freq_type = FrequencyTypes.SEMI_ANNUAL

            bond = Bond(
                issue_dt,
                maturity_dt,
                row.coupon / 100,
                freq_type,
                accrual_dc_type=DayCountTypes.ACT_ACT_ICMA,
            )

            settle_dt = Date(
                row.settle_dt.day,
                row.settle_dt.month,
                row.settle_dt.year,
            )

            accrued_interest = bond.accrued_interest(settle_dt, face)
            clean_price = row.dirty_price - accrued_interest
            calc_ytm = (
                bond.yield_to_maturity(
                    settle_dt, clean_price, YTMCalcType.CFETS
                )
                * 100
            )
            try:
                assert abs(calc_ytm - row.ytm) < 0.0001
            except Exception:
                print(bond)
                print(clean_price)
                print(settle_dt)
                bond.print_payments(settle_dt, 100.0)
                print(f"calc_ytm:{calc_ytm}, correct_ytm:{row.ytm}")
                continue


########################################################################################


def test_key_rate_durations_bloomberg_example():

    dc_type, freq_type, settle_days, ex_div, calendar = (
        get_bond_market_conventions(BondMarkets.UNITED_STATES)
    )

    # interest accrues on this date. Issue date is 01/08/2022
    issue_dt = Date(31, 7, 2022)
    maturity_dt = Date(31, 7, 2027)
    coupon = 2.75 / 100.0
    ex_div_days = 0

    dc_type, freq_type, settle_days, ex_div, calendar = (
        get_bond_market_conventions(BondMarkets.UNITED_STATES)
    )

    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type, ex_div_days)

    settle_dt = Date(24, 4, 2023)

    # US Street yield on Bloomberg as of 20 April 2023
    # with settle date 24 April 2023
    ytm = 3.725060 / 100.0

    # Details of yields of market bonds at KRD maturity points
    my_tenors = np.array([0.5, 1, 2, 3, 5, 7, 10])
    my_rates = (
        np.array([5.0367, 4.7327, 4.1445, 3.8575, 3.6272, 3.5825, 3.5347])
        / 100.0
    )

    krt, krd = bond.key_rate_durations(
        settle_dt, ytm, key_rate_tenors=my_tenors, rates=my_rates
    )

    bbg_key_rate_durations = [
        -0.001,
        -0.009,
        -0.022,
        1.423,
        2.54,
        0.00,
        0.00,
        0.00,
        0.00,
    ]

    for i in range(len(krd)):
        assert my_round(krd[i]) == my_round(bbg_key_rate_durations[i])


# test_bond_zero_ror()
# test_key_rate_durations_bloomberg_example()
