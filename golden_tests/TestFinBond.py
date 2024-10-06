##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import os
import datetime as dt
import pandas as pd
import numpy as np

import sys

sys.path.append("..")

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.calendar import CalendarTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date, from_datetime
from financepy.utils.math import ONE_MILLION
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.bonds.bond_market import get_bond_market_conventions
from financepy.products.bonds.bond_market import BondMarkets
from financepy.products.bonds.bond import YTMCalcType, Bond
from financepy.utils.global_types import SwapTypes


test_cases = FinTestCases(__file__, globalTestCaseMode)


##########################################################################


def build_ibor_curve(value_dt):
    depoDCCType = DayCountTypes.THIRTY_E_360_ISDA
    depos = []
    deposit_rate = 0.050

    depo0 = IborDeposit(value_dt, "1D", deposit_rate, depoDCCType)

    spot_days = 2
    settle_dt = value_dt.add_weekdays(spot_days)

    maturity_dt = settle_dt.add_months(1)
    depo1 = IborDeposit(settle_dt, maturity_dt, deposit_rate, depoDCCType)

    maturity_dt = settle_dt.add_months(3)
    depo2 = IborDeposit(settle_dt, maturity_dt, deposit_rate, depoDCCType)

    maturity_dt = settle_dt.add_months(6)
    depo3 = IborDeposit(settle_dt, maturity_dt, deposit_rate, depoDCCType)

    maturity_dt = settle_dt.add_months(9)
    depo4 = IborDeposit(settle_dt, maturity_dt, deposit_rate, depoDCCType)

    maturity_dt = settle_dt.add_months(12)
    depo5 = IborDeposit(settle_dt, maturity_dt, deposit_rate, depoDCCType)

    depos.append(depo0)
    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)
    depos.append(depo4)
    depos.append(depo5)

    fras = []
    fixed_dcc_type = DayCountTypes.ACT_365F
    fixed_freq_type = FrequencyTypes.SEMI_ANNUAL

    swaps = []

    swap_rate = 0.05
    maturity_dt = settle_dt.add_months(24)
    swap1 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        swap_rate,
        fixed_freq_type,
        fixed_dcc_type,
    )

    #    print(swap1.fixed_leg._payment_dts)

    swaps.append(swap1)

    maturity_dt = settle_dt.add_months(36)
    swap2 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        swap_rate,
        fixed_freq_type,
        fixed_dcc_type,
    )
    swaps.append(swap2)

    #    print(swap2.fixed_leg._payment_dts)

    maturity_dt = settle_dt.add_months(48)
    swap3 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        swap_rate,
        fixed_freq_type,
        fixed_dcc_type,
    )
    swaps.append(swap3)

    #    print(swap3.fixed_leg._payment_dts)

    maturity_dt = settle_dt.add_months(60)
    swap4 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        swap_rate,
        fixed_freq_type,
        fixed_dcc_type,
    )
    swaps.append(swap4)

    #    print(swap4.fixed_leg._payment_dts)

    maturity_dt = settle_dt.add_months(72)
    swap5 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        swap_rate,
        fixed_freq_type,
        fixed_dcc_type,
    )
    swaps.append(swap5)

    #    print(swap5.fixed_leg._payment_dts)

    maturity_dt = settle_dt.add_months(84)
    swap6 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        swap_rate,
        fixed_freq_type,
        fixed_dcc_type,
    )
    swaps.append(swap6)

    #    print(swap6.fixed_leg._payment_dts)

    maturity_dt = settle_dt.add_months(96)
    swap7 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        swap_rate,
        fixed_freq_type,
        fixed_dcc_type,
    )
    swaps.append(swap7)

    #    print(swap7.fixed_leg._payment_dts)

    maturity_dt = settle_dt.add_months(108)
    swap8 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        swap_rate,
        fixed_freq_type,
        fixed_dcc_type,
    )
    swaps.append(swap8)

    #    print(swap8.fixed_leg._payment_dts)

    maturity_dt = settle_dt.add_months(120)
    swap9 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        swap_rate,
        fixed_freq_type,
        fixed_dcc_type,
    )
    swaps.append(swap9)

    #    print(swap9.fixed_leg._payment_dts)

    libor_curve = IborSingleCurve(value_dt, depos, fras, swaps)

    if 1 == 0:
        import numpy as np

        num_steps = 40
        dt = 10 / num_steps
        times = np.linspace(0.0, 10.0, num_steps + 1)

        df0 = 1.0
        for t in times[1:]:
            df1 = libor_curve.df_t(t)
            fwd = (df0 / df1 - 1.0) / dt
            print(t, df1, fwd)
            df0 = df1

    return libor_curve


##########################################################################


def test_bond():

    import pandas as pd

    path = os.path.join(
        os.path.dirname(__file__), ".//data//gilt_bond_prices.txt"
    )
    bond_dataframe = pd.read_csv(path, sep="\t")
    bond_dataframe["mid"] = 0.5 * (
        bond_dataframe["bid"] + bond_dataframe["ask"]
    )

    freq_type = FrequencyTypes.SEMI_ANNUAL
    settle_dt = Date(19, 9, 2012)
    face = ONE_MILLION
    ex_div_days = 0

    for dc_type in DayCountTypes:
        if dc_type == DayCountTypes.ZERO:
            continue
        test_cases.header(
            "MATURITY", "COUPON", "CLEAN_PRICE", "ACCD_DAYS", "ACCRUED", "YTM"
        )

        for _, bond in bond_dataframe.iterrows():
            date_string = bond["maturity"]
            mat_dt_time = dt.datetime.strptime(date_string, "%d-%b-%y")
            maturity_dt = from_datetime(mat_dt_time)
            issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)

            coupon = bond["coupon"] / 100.0
            clean_price = bond["mid"]
            bond = Bond(
                issue_dt, maturity_dt, coupon, freq_type, dc_type, ex_div_days
            )

            ytm = bond.yield_to_maturity(settle_dt, clean_price)
            accrued_int = bond.accrued_int
            accd_days = bond.accrued_days

            test_cases.print(
                "%18s" % maturity_dt,
                "%8.4f" % coupon,
                "%10.4f" % clean_price,
                "%6.0f" % accd_days,
                "%10.4f" % accrued_int,
                "%8.4f" % ytm,
            )

    ###########################################################################
    #  EXAMPLE FROM http://bondtutor.com/btchp4/topic6/topic6.htm

    accrualConvention = DayCountTypes.ACT_ACT_ICMA
    y = 0.062267
    settle_dt = Date(19, 4, 1994)
    issue_dt = Date(15, 7, 1990)
    maturity_dt = Date(15, 7, 1997)
    coupon = 0.085
    ex_div_days = 0
    face = 1000000.0

    freq_type = FrequencyTypes.SEMI_ANNUAL

    bond = Bond(
        issue_dt,
        maturity_dt,
        coupon,
        freq_type,
        accrualConvention,
        ex_div_days,
    )

    test_cases.header("FIELD", "VALUE")
    dirty_price = bond.dirty_price_from_ytm(settle_dt, y)
    test_cases.print("Dirty Price = ", dirty_price)

    clean_price = bond.clean_price_from_ytm(settle_dt, y)
    test_cases.print("Clean Price = ", clean_price)

    accrued_interest = bond.accrued_interest(settle_dt, face)
    test_cases.print("Accrued = ", accrued_interest)

    ytm = bond.yield_to_maturity(settle_dt, clean_price)
    test_cases.print("Yield to Maturity = ", ytm)

    bump = 1e-4
    priceBumpedUp = bond.dirty_price_from_ytm(settle_dt, y + bump)
    test_cases.print("Price Bumped Up:", priceBumpedUp)

    priceBumpedDn = bond.dirty_price_from_ytm(settle_dt, y - bump)
    test_cases.print("Price Bumped Dn:", priceBumpedDn)

    durationByBump = -(priceBumpedUp - dirty_price) / bump
    test_cases.print("Duration by Bump = ", durationByBump)

    duration = bond.dollar_duration(settle_dt, y)
    test_cases.print("Dollar Duration = ", duration)
    test_cases.print("Duration Difference:", duration - durationByBump)

    modified_duration = bond.modified_duration(settle_dt, y)
    test_cases.print("Modified Duration = ", modified_duration)

    macauley_duration = bond.macauley_duration(settle_dt, y)
    test_cases.print("Macauley Duration = ", macauley_duration)

    conv = bond.convexity_from_ytm(settle_dt, y)
    test_cases.print("Convexity = ", conv)

    # ASSET SWAP SPREAD

    # When the libor curve is the flat bond curve then the ASW is zero by
    # definition
    flat_curve = DiscountCurveFlat(settle_dt, ytm, FrequencyTypes.SEMI_ANNUAL)

    test_cases.header("FIELD", "VALUE")

    clean_price = bond.clean_price_from_ytm(settle_dt, ytm)
    asw = bond.asset_swap_spread(settle_dt, clean_price, flat_curve)
    test_cases.print("Discounted on Bond Curve ASW:", asw * 10000)

    # When the libor curve is the Libor curve then the ASW is positive
    libor_curve = build_ibor_curve(settle_dt)
    asw = bond.asset_swap_spread(settle_dt, clean_price, libor_curve)
    oas = bond.option_adjusted_spread(settle_dt, clean_price, libor_curve)
    test_cases.print("Discounted on LIBOR Curve ASW:", asw * 10000)
    test_cases.print("Discounted on LIBOR Curve OAS:", oas * 10000)

    p = 90.0
    asw = bond.asset_swap_spread(settle_dt, p, libor_curve)
    oas = bond.option_adjusted_spread(settle_dt, p, libor_curve)
    test_cases.print("Deep discount bond at 90 ASW:", asw * 10000)
    test_cases.print("Deep discount bond at 90 OAS:", oas * 10000)

    p = 100.0
    asw = bond.asset_swap_spread(settle_dt, p, libor_curve)
    oas = bond.option_adjusted_spread(settle_dt, p, libor_curve)
    test_cases.print("Par bond at 100 ASW:", asw * 10000)
    test_cases.print("Par bond at 100 OAS:", oas * 10000)

    p = 120.0
    asw = bond.asset_swap_spread(settle_dt, p, libor_curve)
    oas = bond.option_adjusted_spread(settle_dt, p, libor_curve)
    test_cases.print("Above par bond at 120 ASW:", asw * 10000)
    test_cases.print("Above par bond at 120 OAS:", oas * 10000)

    ##########################################################################
    # https://data.bloomberglp.com/bat/sites/3/2017/07/SF-2017_Paul-Fjeldsted.pdf
    # Page 10 TREASURY NOTE SCREENSHOT
    ##########################################################################

    test_cases.banner("BLOOMBERG US TREASURY EXAMPLE")
    settle_dt = Date(21, 7, 2017)
    issue_dt = Date(15, 5, 2010)
    maturity_dt = Date(15, 5, 2027)
    coupon = 0.02375
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    ex_div_days = 0
    face = 1000000.0

    bond = Bond(
        issue_dt,
        maturity_dt,
        coupon,
        freq_type,
        dc_type,
        ex_div_days,
        CalendarTypes.UNITED_STATES,
    )

    test_cases.header("FIELD", "VALUE")
    clean_price = 99.7808417

    yld = bond.current_yield(clean_price)
    test_cases.print("Current Yield = ", yld)

    ytm = bond.yield_to_maturity(settle_dt, clean_price, YTMCalcType.UK_DMO)
    test_cases.print("UK DMO Yield To Maturity = ", ytm)

    ytm = bond.yield_to_maturity(settle_dt, clean_price, YTMCalcType.US_STREET)
    test_cases.print("US STREET Yield To Maturity = ", ytm)

    ytm = bond.yield_to_maturity(
        settle_dt, clean_price, YTMCalcType.US_TREASURY
    )
    test_cases.print("US TREASURY Yield To Maturity = ", ytm)

    dirty_price = bond.dirty_price_from_ytm(
        settle_dt, ytm, YTMCalcType.US_TREASURY
    )

    test_cases.print("Dirty Price = ", dirty_price)

    clean_price = bond.clean_price_from_ytm(
        settle_dt, ytm, YTMCalcType.US_TREASURY
    )
    test_cases.print("Clean Price = ", clean_price)

    accrued_interest = bond.accrued_interest(settle_dt, face)
    test_cases.print("Accrued = ", accrued_interest)

    accddays = bond.accrued_days
    test_cases.print("Accrued Days = ", accddays)

    duration = bond.dollar_duration(settle_dt, ytm, YTMCalcType.US_STREET)
    test_cases.print("Dollar Duration = ", duration)

    modified_duration = bond.modified_duration(settle_dt, ytm)
    test_cases.print("Modified Duration = ", modified_duration)

    macauley_duration = bond.macauley_duration(settle_dt, ytm)
    test_cases.print("Macauley Duration = ", macauley_duration)

    conv = bond.convexity_from_ytm(settle_dt, ytm)
    test_cases.print("Convexity = ", conv)

    ##########################################################################
    # Page 11 APPLE NOTE SCREENSHOT
    ##########################################################################

    test_cases.banner("BLOOMBERG APPLE CORP BOND EXAMPLE")
    settle_dt = Date(21, 7, 2017)
    issue_dt = Date(13, 5, 2012)
    maturity_dt = Date(13, 5, 2022)
    coupon = 0.027
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.THIRTY_E_360_ISDA
    ex_div_days = 0
    face = 100.0

    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type, ex_div_days)

    test_cases.header("FIELD", "VALUE")
    clean_price = 101.581564

    yld = bond.current_yield(clean_price)
    test_cases.print("Current Yield", yld)

    ytm = bond.yield_to_maturity(settle_dt, clean_price, YTMCalcType.UK_DMO)
    test_cases.print("UK DMO Yield To Maturity", ytm)

    ytm = bond.yield_to_maturity(settle_dt, clean_price, YTMCalcType.US_STREET)
    test_cases.print("US STREET Yield To Maturity", ytm)

    ytm = bond.yield_to_maturity(
        settle_dt, clean_price, YTMCalcType.US_TREASURY
    )
    test_cases.print("US TREASURY Yield To Maturity", ytm)

    dirty_price = bond.dirty_price_from_ytm(settle_dt, ytm)
    test_cases.print("Dirty Price", dirty_price)

    clean_price = bond.clean_price_from_ytm(settle_dt, ytm)
    test_cases.print("Clean Price", clean_price)

    accddays = bond.accrued_days
    test_cases.print("Accrued Days", accddays)

    accrued_interest = bond.accrued_interest(settle_dt, face)
    test_cases.print("Accrued", accrued_interest)

    duration = bond.dollar_duration(settle_dt, ytm)
    test_cases.print("Dollar Duration", duration)

    modified_duration = bond.modified_duration(settle_dt, ytm)
    test_cases.print("Modified Duration", modified_duration)

    macauley_duration = bond.macauley_duration(settle_dt, ytm)
    test_cases.print("Macauley Duration", macauley_duration)

    conv = bond.convexity_from_ytm(settle_dt, ytm)
    test_cases.print("Convexity", conv)


###############################################################################


def test_bond_ex_dividend():

    issue_dt = Date(7, 9, 2000)
    maturity_dt = Date(7, 9, 2020)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    face = 100.0
    ex_div_days = 7
    test_cases.header("LABEL", "VALUE")

    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type, ex_div_days)
    settle_dt = Date(7, 9, 2003)
    accrued = bond.accrued_interest(settle_dt, face)

    test_cases.print("settle_dt:", settle_dt)
    test_cases.print("Accrued:", accrued)

    ###########################################################################
    test_cases.banner(
        "======================================================="
    )
    test_cases.header("SETTLEMENT", "DIRTY PRICE", "ACCRUED", "CLEAN PRICE")

    issue_dt = Date(7, 9, 2000)
    maturity_dt = Date(7, 9, 2020)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    face = 100.0
    ex_div_days = 7

    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type, ex_div_days)

    settle_dt = Date(25, 8, 2010)

    ytm = 0.05

    for _ in range(0, 13):
        settle_dt = settle_dt.add_days(1)
        accrued = bond.accrued_interest(settle_dt, face)
        dirty_price = bond.dirty_price_from_ytm(settle_dt, ytm)
        clean_price = dirty_price - accrued
        test_cases.print(settle_dt, dirty_price, accrued, clean_price)


###############################################################################


def test_bond_payment_dates():

    from financepy.products.bonds.bond import Bond
    from financepy.utils import Date, DayCountTypes, FrequencyTypes

    bond = Bond(
        issue_dt=Date(7, 6, 2021),
        maturity_dt=Date(7, 6, 2031),
        coupon=0.0341,
        freq_type=FrequencyTypes.ANNUAL,
        dc_type=DayCountTypes.ACT_ACT_ISDA,
    )
    bond._calculate_payment_dts()

    if 1 == 0:
        print(bond.flow_amounts)
        print(bond.cpn_dts)
        print(bond._payment_dts)


###############################################################################


def test_bond_ror():

    path = os.path.join(
        os.path.dirname(__file__), ".//data//test_cases_bond_ror.csv"
    )
    df = pd.read_csv(path, parse_dates=["buy_date", "sell_date"])
    # A 10-year bond with 1 coupon per year. code: 210215

    bond = Bond(
        issue_dt=Date(13, 9, 2021),
        maturity_dt=Date(13, 9, 2031),
        coupon=0.0312,
        freq_type=FrequencyTypes.ANNUAL,
        dc_type=DayCountTypes.ACT_ACT_ICMA,
    )

    test_cases.header(
        "bond_code",
        "buy_date",
        "buy_ytm",
        "buy_price",
        "sell_date",
        "sell_ytm",
        "sell_price",
        "simple_return",
        "irr",
    )

    for row in df.itertuples(index=False):

        buy_date = Date(
            row.buy_date.day, row.buy_date.month, row.buy_date.year
        )
        sell_date = Date(
            row.sell_date.day, row.sell_date.month, row.sell_date.year
        )
        buy_price = bond.dirty_price_from_ytm(
            buy_date, row.buy_ytm, YTMCalcType.US_STREET
        )
        sell_price = bond.dirty_price_from_ytm(
            sell_date, row.sell_ytm, YTMCalcType.US_STREET
        )
        simple, irr, pnl = bond.calc_ror(
            buy_date, sell_date, row.buy_ytm, row.sell_ytm
        )

        test_cases.print(
            row.bond_code,
            buy_date,
            row.buy_ytm,
            buy_price,
            sell_date,
            row.sell_ytm,
            sell_price,
            simple,
            irr,
        )


###############################################################################


def test_bond_eom():

    # Bonds that mature on an EOM date have flows on EOM dates
    issue_dt = Date(30, 11, 2022)
    settle_dt = Date(6, 2, 2023)
    maturity_dt = Date(30, 11, 2024)
    coupon = 0.045
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    ex_div_days = 0

    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type, ex_div_days)

    ai = bond.accrued_interest(settle_dt)  # should be 8406.593406


###############################################################################


def test_key_rate_durations():

    issue_dt = Date(31, 7, 2022)
    maturity_dt = Date(31, 7, 2027)
    coupon = 0.0275
    ex_div_days = 0

    dc_type, freq_type, settle_days, exDiv, calendar = (
        get_bond_market_conventions(BondMarkets.UNITED_STATES)
    )

    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type, ex_div_days)

    settle_dt = Date(24, 4, 2023)

    ytm = 3.725060 / 100.0

    krt, krd = bond.key_rate_durations(settle_dt, ytm)


#    print(key_rate_tenors)
#    print(key_rate_durations)

###############################################################################


def test_key_rate_durations_bloomberg_example():

    dc_type, frequencyType, settle_days, exDiv, calendar = (
        get_bond_market_conventions(BondMarkets.UNITED_STATES)
    )

    # interest accrues on this date. Issue date is 01/08/2022
    issue_dt = Date(31, 7, 2022)
    maturity_dt = Date(31, 7, 2027)
    coupon = 2.75 / 100.0
    ex_div_days = 0

    dc_type, freq_type, settle_days, exDiv, calendar = (
        get_bond_market_conventions(BondMarkets.UNITED_STATES)
    )

    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type, ex_div_days)

    settle_dt = Date(24, 4, 2023)

    # US Street yield on Bloomberg as of 20 April 2023
    # with settle date 24 April 2023
    ytm = 3.725060 / 100

    # Details of yields of market bonds at KRD maturity points
    my_tenors = np.array([0.5, 1, 2, 3, 5, 7, 10])

    my_rates = (
        np.array([5.0367, 4.7327, 4.1445, 3.8575, 3.6272, 3.5825, 3.5347])
        / 100
    )

    krt, krd = bond.key_rate_durations(
        settle_dt, ytm, key_rate_tenors=my_tenors, rates=my_rates
    )


#    print(key_rate_tenors)
#    print(key_rate_durations)

# Differences due to bonds not sitting exactly on these maturity points ?
# Did BBG interpolate ?

###############################################################################


def test_oas():

    issue_dt = Date(15, 5, 2010)
    maturity_dt = Date(15, 5, 2027)
    coupon = 0.02375
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA

    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    liborFlatRate = 0.0275
    settle_dt = Date(21, 7, 2017)

    liborFlatCurve = DiscountCurveFlat(
        settle_dt, liborFlatRate, FrequencyTypes.SEMI_ANNUAL
    )

    # I specified face to be 100 - if face is 1 then this must be 0.99780842
    clean_price = 99.780842

    oas = (
        bond.option_adjusted_spread(settle_dt, clean_price, liborFlatCurve)
        * 10000
    )

    if (oas - (-34.95)) > 0.01:
        print("OAS incorrect")


###############################################################################


def test_div_dts():

    issue_dt = Date(15, 5, 2020)
    maturity_dt = Date(15, 5, 2035)
    coupon = 0.02375
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    face = 125000
    ex_div_days = 10

    bond = Bond(
        issue_dt, maturity_dt, coupon, freq_type, accrual_type, ex_div_days
    )

    print(bond)

    clean_price = 99.7808417  # if face is 1 then this must be 0.99780842

    settle_dt = Date(15, 5, 2023)
    print(bond.bond_payments(settle_dt, face))

    current_yield = bond.current_yield(clean_price) * 100
    print("Currnt Yield: %10.5f %%" % (current_yield))

    ytm = bond.yield_to_maturity(settle_dt, clean_price) * 100.0
    print("Yield to Mat: %10.5f %%" % (ytm))


###############################################################################

test_bond()
test_oas()
test_bond_ex_dividend()
test_bond_payment_dates()
test_bond_ror()
test_bond_eom()
test_key_rate_durations()
test_key_rate_durations_bloomberg_example()

test_cases.compareTestCases()
