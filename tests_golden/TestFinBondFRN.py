###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.global_types import SwapTypes
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_swap import IborSwap
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.bonds.bond_frn import BondFRN
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.utils.date import Date
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def build_Ibor_Curve(valuation_date):

    depoDCCType = DayCountTypes.THIRTY_E_360_ISDA
    depos = []

    payFixed = SwapTypes.PAY

    spot_days = 2
    settlement_date = valuation_date.add_weekdays(spot_days)

    deposit_rate = 0.050
    maturity_date = settlement_date.add_months(1)
    depo1 = IborDeposit(
        settlement_date,
        maturity_date,
        deposit_rate,
        depoDCCType)

    maturity_date = settlement_date.add_months(3)
    depo2 = IborDeposit(
        settlement_date,
        maturity_date,
        deposit_rate,
        depoDCCType)

    maturity_date = settlement_date.add_months(6)
    depo3 = IborDeposit(
        settlement_date,
        maturity_date,
        deposit_rate,
        depoDCCType)

    maturity_date = settlement_date.add_months(9)
    depo4 = IborDeposit(
        settlement_date,
        maturity_date,
        deposit_rate,
        depoDCCType)

    maturity_date = settlement_date.add_months(12)
    depo5 = IborDeposit(
        settlement_date,
        maturity_date,
        deposit_rate,
        depoDCCType)

    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)
    depos.append(depo4)
    depos.append(depo5)

    fras = []
    fixedDCCType = DayCountTypes.ACT_365F
    fixedFreqType = FrequencyTypes.SEMI_ANNUAL

    swaps = []

    swap_rate = 0.05
    maturity_date = settlement_date.add_months(24)
    swap1 = IborSwap(
        settlement_date,
        maturity_date,
        swap_rate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap1)

    maturity_date = settlement_date.add_months(36)
    swap2 = IborSwap(
        settlement_date,
        maturity_date,
        swap_rate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap2)

    maturity_date = settlement_date.add_months(48)
    swap3 = IborSwap(
        settlement_date,
        maturity_date,
        swap_rate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap3)

    maturity_date = settlement_date.add_months(60)
    swap4 = IborSwap(
        settlement_date,
        maturity_date,
        swap_rate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap4)

    maturity_date = settlement_date.add_months(72)
    swap5 = IborSwap(
        settlement_date,
        maturity_date,
        swap_rate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap5)

    maturity_date = settlement_date.add_months(84)
    swap6 = IborSwap(
        settlement_date,
        maturity_date,
        swap_rate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap6)

    maturity_date = settlement_date.add_months(96)
    swap7 = IborSwap(
        settlement_date,
        maturity_date,
        swap_rate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap7)

    maturity_date = settlement_date.add_months(108)
    swap8 = IborSwap(
        settlement_date,
        maturity_date,
        swap_rate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap8)

    maturity_date = settlement_date.add_months(120)
    swap9 = IborSwap(
        settlement_date,
        maturity_date,
        swap_rate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap9)

    libor_curve = IborSingleCurve(valuation_date,
                                  depos,
                                  fras,
                                  swaps)

    if 1 == 0:
        import numpy as np
        num_steps = 40
        dt = 10 / num_steps
        times = np.linspace(0.0, 10.0, num_steps + 1)

        df0 = 1.0
        for t in times[1:]:
            df1 = libor_curve.df(t)
            fwd = (df0 / df1 - 1.0) / dt
            print(t, df1, fwd)
            df0 = df1

    return libor_curve

##########################################################################


def test_BondFRN():

    # https://data.bloomberglp.com/bat/sites/3/2017/07/SF-2017_Paul-Fjeldsted.pdf
    # I have a day out problem on the accrued interest - should be 71 and not 72 days
    # Other than that agreement on the DM is very good.

    ##########################################################################
    # CITIGROUP FRN SCREENSHOT
    ##########################################################################

    testCases.banner("BLOOMBERG CITIGROUP FRN EXAMPLE")
    issue_date = Date(10, 11, 2010)
    maturity_date = Date(10, 11, 2021)
    quoted_margin = 0.0025
    freq_type = FrequencyTypes.QUARTERLY
    accrual_type = DayCountTypes.THIRTY_E_360
    face = 1000000

    bond = BondFRN(issue_date,
                   maturity_date,
                   quoted_margin,
                   freq_type,
                   accrual_type,
                   face)

    testCases.header("FIELD", "VALUE")
    clean_price = 96.793
    resetIbor = 0.0143456 - quoted_margin
    current_ibor = 0.0120534
    future_ibors = 0.0130522

    settlement_date = Date(21, 7, 2017)

    dm = bond.discount_margin(settlement_date,
                              resetIbor,
                              current_ibor,
                              future_ibors,
                              clean_price)

    testCases.print("Discount Margin (bp) = ", dm * 10000)

    full_price = bond.full_price_from_dm(settlement_date,
                                         resetIbor,
                                         current_ibor,
                                         future_ibors,
                                         dm)

    testCases.print("Full Price = ", full_price)

    lastCouponDt = bond._pcd
    testCases.print("Last Coupon Date = ", str(lastCouponDt))

    accddays = bond._accrued_days
    testCases.print("Accrued Days = ", accddays)

    accdAmount = bond._accrued_interest
    testCases.print("Accrued Amount = ", accdAmount)

    principal = bond.principal(settlement_date,
                               resetIbor,
                               current_ibor,
                               future_ibors,
                               dm)

    testCases.print("Dollar Principal = ", principal)

    duration = bond.dollar_duration(settlement_date,
                                    resetIbor,
                                    current_ibor,
                                    future_ibors,
                                    dm)

    testCases.print("Dollar Rate Duration = ", duration)

    modified_duration = bond.modified_duration(settlement_date,
                                               resetIbor,
                                               current_ibor,
                                               future_ibors,
                                               dm)

    testCases.print("Modified Rate Duration = ", modified_duration)

    macauley_duration = bond.macauley_duration(settlement_date,
                                               resetIbor,
                                               current_ibor,
                                               future_ibors,
                                               dm)

    testCases.print("Macauley Duration = ", macauley_duration)

    convexity = bond.convexity_from_dm(settlement_date,
                                       resetIbor,
                                       current_ibor,
                                       future_ibors,
                                       dm)

    testCases.print("Convexity = ", convexity)

    duration = bond.dollar_credit_duration(settlement_date,
                                           resetIbor,
                                           current_ibor,
                                           future_ibors,
                                           dm)

    testCases.print("Dollar Credit Duration = ", duration)

    modified_duration = bond.modified_credit_duration(settlement_date,
                                                      resetIbor,
                                                      current_ibor,
                                                      future_ibors,
                                                      dm)

    testCases.print("Modified Credit Duration = ", modified_duration)

##########################################################################
# EXAMPLE
# https://ebrary.net/14293/economics/actual_floater
##########################################################################

    testCases.banner("BLOOMBERG CITIGROUP FRN EXAMPLE II")
    issue_date = Date(28, 3, 2000)
    settlement_date = Date(28, 3, 2014)
    maturity_date = Date(3, 2, 2021)
    quoted_margin = 0.0020
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.THIRTY_E_360_ISDA
    face = 1000000.0

    bond = BondFRN(issue_date,
                   maturity_date,
                   quoted_margin,
                   freq_type,
                   accrual_type,
                   face)

    testCases.header("FIELD", "VALUE")
    clean_price = 93.08
    resetIbor = 0.00537 - quoted_margin
    current_ibor = 0.027558
    future_ibors = 0.03295

    dm = bond.discount_margin(settlement_date,
                              resetIbor,
                              current_ibor,
                              future_ibors,
                              clean_price)

    testCases.print("Discount Margin (bp) = ", dm * 10000)

    full_price = bond.full_price_from_dm(settlement_date,
                                         resetIbor,
                                         current_ibor,
                                         future_ibors,
                                         dm)

    testCases.print("Full Price = ", full_price)

    lastCouponDt = bond._pcd
    testCases.print("Last Coupon Date = ", str(lastCouponDt))

    accddays = bond._accrued_days
    testCases.print("Accrued Days = ", accddays)

    accdAmount = bond._accrued_interest
    testCases.print("Accrued Amount = ", accdAmount)

    principal = bond.principal(settlement_date,
                               resetIbor,
                               current_ibor,
                               future_ibors,
                               dm)

    testCases.print("Dollar Principal = ", principal)

    duration = bond.dollar_duration(settlement_date,
                                    resetIbor,
                                    current_ibor,
                                    future_ibors,
                                    dm)

    testCases.print("Dollar Rate Duration = ", duration)

    modified_duration = bond.modified_duration(settlement_date,
                                               resetIbor,
                                               current_ibor,
                                               future_ibors,
                                               dm)

    testCases.print("Modified Rate Duration = ", modified_duration)

    macauley_duration = bond.macauley_duration(settlement_date,
                                               resetIbor,
                                               current_ibor,
                                               future_ibors,
                                               dm)

    testCases.print("Macauley Duration = ", macauley_duration)

    convexity = bond.convexity_from_dm(settlement_date,
                                       resetIbor,
                                       current_ibor,
                                       future_ibors,
                                       dm)

    testCases.print("Convexity = ", convexity)

    principal = bond.principal(settlement_date,
                               resetIbor,
                               current_ibor,
                               future_ibors,
                               dm)

    testCases.print("Principal = ", principal)

    duration = bond.dollar_credit_duration(settlement_date,
                                           resetIbor,
                                           current_ibor,
                                           future_ibors,
                                           dm)

    testCases.print("Dollar Credit Duration = ", duration)

    modified_duration = bond.modified_credit_duration(settlement_date,
                                                      resetIbor,
                                                      current_ibor,
                                                      future_ibors,
                                                      dm)

    testCases.print("Modified Credit Duration = ", modified_duration)

##########################################################################


test_BondFRN()
testCases.compareTestCases()
