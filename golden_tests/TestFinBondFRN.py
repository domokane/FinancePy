###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys

sys.path.append("..")

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.global_types import SwapTypes
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_swap import IborSwap
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.bonds.bond_frn import BondFRN
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.utils.date import Date


test_cases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def build_Ibor_Curve(value_dt):

    depoDCCType = DayCountTypes.THIRTY_E_360_ISDA
    depos = []

    pay_fixed = SwapTypes.PAY

    spot_days = 2
    settle_dt = value_dt.add_weekdays(spot_days)

    deposit_rate = 0.050
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
        swap_rate,
        pay_fixed,
        fixed_freq_type,
        fixed_dcc_type,
    )
    swaps.append(swap1)

    maturity_dt = settle_dt.add_months(36)
    swap2 = IborSwap(
        settle_dt,
        maturity_dt,
        swap_rate,
        pay_fixed,
        fixed_freq_type,
        fixed_dcc_type,
    )
    swaps.append(swap2)

    maturity_dt = settle_dt.add_months(48)
    swap3 = IborSwap(
        settle_dt,
        maturity_dt,
        swap_rate,
        pay_fixed,
        fixed_freq_type,
        fixed_dcc_type,
    )
    swaps.append(swap3)

    maturity_dt = settle_dt.add_months(60)
    swap4 = IborSwap(
        settle_dt,
        maturity_dt,
        swap_rate,
        pay_fixed,
        fixed_freq_type,
        fixed_dcc_type,
    )
    swaps.append(swap4)

    maturity_dt = settle_dt.add_months(72)
    swap5 = IborSwap(
        settle_dt,
        maturity_dt,
        swap_rate,
        pay_fixed,
        fixed_freq_type,
        fixed_dcc_type,
    )
    swaps.append(swap5)

    maturity_dt = settle_dt.add_months(84)
    swap6 = IborSwap(
        settle_dt,
        maturity_dt,
        swap_rate,
        pay_fixed,
        fixed_freq_type,
        fixed_dcc_type,
    )
    swaps.append(swap6)

    maturity_dt = settle_dt.add_months(96)
    swap7 = IborSwap(
        settle_dt,
        maturity_dt,
        swap_rate,
        pay_fixed,
        fixed_freq_type,
        fixed_dcc_type,
    )
    swaps.append(swap7)

    maturity_dt = settle_dt.add_months(108)
    swap8 = IborSwap(
        settle_dt,
        maturity_dt,
        swap_rate,
        pay_fixed,
        fixed_freq_type,
        fixed_dcc_type,
    )
    swaps.append(swap8)

    maturity_dt = settle_dt.add_months(120)
    swap9 = IborSwap(
        settle_dt,
        maturity_dt,
        swap_rate,
        pay_fixed,
        fixed_freq_type,
        fixed_dcc_type,
    )
    swaps.append(swap9)

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


def test_BondFRN():

    # https://data.bloomberglp.com/bat/sites/3/2017/07/SF-2017_Paul-Fjeldsted.pdf
    # I have a day out problem on the accrued interest - should be 71 and
    # not 72 days
    # Other than that agreement on the DM is very good.

    ##########################################################################
    # CITIGROUP FRN SCREENSHOT
    ##########################################################################

    test_cases.banner("BLOOMBERG CITIGROUP FRN EXAMPLE")
    issue_dt = Date(10, 11, 2010)
    maturity_dt = Date(10, 11, 2021)
    quoted_margin = 0.0025
    freq_type = FrequencyTypes.QUARTERLY
    dc_type = DayCountTypes.THIRTY_E_360

    bond = BondFRN(issue_dt, maturity_dt, quoted_margin, freq_type, dc_type)

    test_cases.header("FIELD", "VALUE")
    clean_price = 96.793
    resetIbor = 0.0143456 - quoted_margin
    current_ibor = 0.0120534
    future_ibors = 0.0130522

    settle_dt = Date(21, 7, 2017)

    dm = bond.discount_margin(
        settle_dt, resetIbor, current_ibor, future_ibors, clean_price
    )

    test_cases.print("Discount Margin (bp) = ", dm * 10000)

    dirty_price = bond.dirty_price_from_dm(
        settle_dt, resetIbor, current_ibor, future_ibors, dm
    )

    test_cases.print("Dirty Price = ", dirty_price)

    lastCouponDt = bond.pcd
    test_cases.print("Last Coupon Date = ", str(lastCouponDt))

    accddays = bond.accrued_days
    test_cases.print("Accrued Days = ", accddays)

    accdAmount = bond.accrued_int
    test_cases.print("Accrued Amount = ", accdAmount)

    principal = bond.principal(
        settle_dt, resetIbor, current_ibor, future_ibors, dm
    )

    test_cases.print("Dollar Principal = ", principal)

    duration = bond.dollar_duration(
        settle_dt, resetIbor, current_ibor, future_ibors, dm
    )

    test_cases.print("Dollar Rate Duration = ", duration)

    modified_duration = bond.modified_duration(
        settle_dt, resetIbor, current_ibor, future_ibors, dm
    )

    test_cases.print("Modified Rate Duration = ", modified_duration)

    macauley_duration = bond.macauley_duration(
        settle_dt, resetIbor, current_ibor, future_ibors, dm
    )

    test_cases.print("Macauley Duration = ", macauley_duration)

    convexity = bond.convexity_from_dm(
        settle_dt, resetIbor, current_ibor, future_ibors, dm
    )

    test_cases.print("Convexity = ", convexity)

    duration = bond.dollar_credit_duration(
        settle_dt, resetIbor, current_ibor, future_ibors, dm
    )

    test_cases.print("Dollar Credit Duration = ", duration)

    modified_duration = bond.modified_credit_duration(
        settle_dt, resetIbor, current_ibor, future_ibors, dm
    )

    test_cases.print("Modified Credit Duration = ", modified_duration)

    ##########################################################################
    # EXAMPLE
    # https://ebrary.net/14293/economics/actual_floater
    ##########################################################################

    test_cases.banner("BLOOMBERG CITIGROUP FRN EXAMPLE II")
    issue_dt = Date(28, 3, 2000)
    settle_dt = Date(28, 3, 2014)
    maturity_dt = Date(3, 2, 2021)
    quoted_margin = 0.0020
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.THIRTY_E_360_ISDA

    bond = BondFRN(issue_dt, maturity_dt, quoted_margin, freq_type, dc_type)

    test_cases.header("FIELD", "VALUE")
    clean_price = 93.08
    resetIbor = 0.00537 - quoted_margin
    current_ibor = 0.027558
    future_ibors = 0.03295

    dm = bond.discount_margin(
        settle_dt, resetIbor, current_ibor, future_ibors, clean_price
    )

    test_cases.print("Discount Margin (bp) = ", dm * 10000)

    dirty_price = bond.dirty_price_from_dm(
        settle_dt, resetIbor, current_ibor, future_ibors, dm
    )

    test_cases.print("Dirty Price = ", dirty_price)

    lastCouponDt = bond.pcd
    test_cases.print("Last Coupon Date = ", str(lastCouponDt))

    accddays = bond.accrued_days
    test_cases.print("Accrued Days = ", accddays)

    accdAmount = bond.accrued_int
    test_cases.print("Accrued Amount = ", accdAmount)

    principal = bond.principal(
        settle_dt, resetIbor, current_ibor, future_ibors, dm
    )

    test_cases.print("Dollar Principal = ", principal)

    duration = bond.dollar_duration(
        settle_dt, resetIbor, current_ibor, future_ibors, dm
    )

    test_cases.print("Dollar Rate Duration = ", duration)

    modified_duration = bond.modified_duration(
        settle_dt, resetIbor, current_ibor, future_ibors, dm
    )

    test_cases.print("Modified Rate Duration = ", modified_duration)

    macauley_duration = bond.macauley_duration(
        settle_dt, resetIbor, current_ibor, future_ibors, dm
    )

    test_cases.print("Macauley Duration = ", macauley_duration)

    convexity = bond.convexity_from_dm(
        settle_dt, resetIbor, current_ibor, future_ibors, dm
    )

    test_cases.print("Convexity = ", convexity)

    principal = bond.principal(
        settle_dt, resetIbor, current_ibor, future_ibors, dm
    )

    test_cases.print("Principal = ", principal)

    duration = bond.dollar_credit_duration(
        settle_dt, resetIbor, current_ibor, future_ibors, dm
    )

    test_cases.print("Dollar Credit Duration = ", duration)

    modified_duration = bond.modified_credit_duration(
        settle_dt, resetIbor, current_ibor, future_ibors, dm
    )

    test_cases.print("Modified Credit Duration = ", modified_duration)


##########################################################################


test_BondFRN()
test_cases.compareTestCases()
