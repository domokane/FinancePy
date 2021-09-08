##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

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
from financepy.products.bonds.bond import Bond
from financepy.products.bonds.bond import YTMCalcType
from financepy.utils.global_types import SwapTypes
import os
import datetime as dt

import sys

sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)


##########################################################################

def build_Ibor_Curve(valuation_date):
    depoDCCType = DayCountTypes.THIRTY_E_360_ISDA
    depos = []
    deposit_rate = 0.050

    depo0 = IborDeposit(
        valuation_date,
        "1D",
        deposit_rate,
        depoDCCType)

    spot_days = 2
    settlement_date = valuation_date.add_weekdays(spot_days)

    maturity_date = settlement_date.add_months(1)
    depo1 = IborDeposit(settlement_date,
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

    depos.append(depo0)
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
        SwapTypes.PAY,
        swap_rate,
        fixedFreqType,
        fixedDCCType)

    #    print(swap1._fixed_leg._payment_dates)

    swaps.append(swap1)

    maturity_date = settlement_date.add_months(36)
    swap2 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        swap_rate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap2)

    #    print(swap2._fixed_leg._payment_dates)

    maturity_date = settlement_date.add_months(48)
    swap3 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        swap_rate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap3)

    #    print(swap3._fixed_leg._payment_dates)

    maturity_date = settlement_date.add_months(60)
    swap4 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        swap_rate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap4)

    #    print(swap4._fixed_leg._payment_dates)

    maturity_date = settlement_date.add_months(72)
    swap5 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        swap_rate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap5)

    #    print(swap5._fixed_leg._payment_dates)

    maturity_date = settlement_date.add_months(84)
    swap6 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        swap_rate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap6)

    #    print(swap6._fixed_leg._payment_dates)

    maturity_date = settlement_date.add_months(96)
    swap7 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        swap_rate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap7)

    #    print(swap7._fixed_leg._payment_dates)

    maturity_date = settlement_date.add_months(108)
    swap8 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        swap_rate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap8)

    #    print(swap8._fixed_leg._payment_dates)

    maturity_date = settlement_date.add_months(120)
    swap9 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        swap_rate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap9)

    #    print(swap9._fixed_leg._payment_dates)

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


def test_Bond():
    import pandas as pd
    path = os.path.join(os.path.dirname(__file__), './data/giltBondPrices.txt')
    bondDataFrame = pd.read_csv(path, sep='\t')
    bondDataFrame['mid'] = 0.5 * (bondDataFrame['bid'] + bondDataFrame['ask'])

    freq_type = FrequencyTypes.SEMI_ANNUAL
    settlement_date = Date(19, 9, 2012)
    face = ONE_MILLION

    for accrual_type in DayCountTypes:

        testCases.header("MATURITY", "COUPON", "CLEAN_PRICE", "ACCD_DAYS",
                         "ACCRUED", "YTM")

        for _, bond in bondDataFrame.iterrows():
            date_string = bond['maturity']
            matDatetime = dt.datetime.strptime(date_string, '%d-%b-%y')
            maturityDt = from_datetime(matDatetime)
            issueDt = Date(maturityDt._d, maturityDt._m, 2000)

            coupon = bond['coupon'] / 100.0
            clean_price = bond['mid']
            bond = Bond(issueDt, maturityDt,
                        coupon, freq_type, accrual_type, 100)

            ytm = bond.yield_to_maturity(settlement_date, clean_price)
            accrued_interest = bond._accrued_interest
            accd_days = bond._accrued_days

            testCases.print("%18s" % maturityDt, "%8.4f" % coupon,
                            "%10.4f" % clean_price, "%6.0f" % accd_days,
                            "%10.4f" % accrued_interest, "%8.4f" % ytm)

    ###########################################################################
    #  EXAMPLE FROM http://bondtutor.com/btchp4/topic6/topic6.htm

    accrualConvention = DayCountTypes.ACT_ACT_ICMA
    y = 0.062267
    settlement_date = Date(19, 4, 1994)
    issue_date = Date(15, 7, 1990)
    maturity_date = Date(15, 7, 1997)
    coupon = 0.085
    face = ONE_MILLION
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issue_date, maturity_date,
                coupon, freq_type, accrualConvention, face)

    testCases.header("FIELD", "VALUE")
    full_price = bond.full_price_from_ytm(settlement_date, y)
    testCases.print("Full Price = ", full_price)
    clean_price = bond.clean_price_from_ytm(settlement_date, y)
    testCases.print("Clean Price = ", clean_price)
    accrued_interest = bond._accrued_interest
    testCases.print("Accrued = ", accrued_interest)
    ytm = bond.yield_to_maturity(settlement_date, clean_price)
    testCases.print("Yield to Maturity = ", ytm)

    bump = 1e-4
    priceBumpedUp = bond.full_price_from_ytm(settlement_date, y + bump)
    testCases.print("Price Bumped Up:", priceBumpedUp)

    priceBumpedDn = bond.full_price_from_ytm(settlement_date, y - bump)
    testCases.print("Price Bumped Dn:", priceBumpedDn)

    durationByBump = -(priceBumpedUp - full_price) / bump
    testCases.print("Duration by Bump = ", durationByBump)

    duration = bond.dollar_duration(settlement_date, y)
    testCases.print("Dollar Duration = ", duration)
    testCases.print("Duration Difference:", duration - durationByBump)

    modified_duration = bond.modified_duration(settlement_date, y)
    testCases.print("Modified Duration = ", modified_duration)

    macauley_duration = bond.macauley_duration(settlement_date, y)
    testCases.print("Macauley Duration = ", macauley_duration)

    conv = bond.convexity_from_ytm(settlement_date, y)
    testCases.print("Convexity = ", conv)

    # ASSET SWAP SPREAD

    # When the libor curve is the flat bond curve then the ASW is zero by
    # definition
    flat_curve = DiscountCurveFlat(settlement_date,
                                   ytm,
                                   FrequencyTypes.SEMI_ANNUAL)

    testCases.header("FIELD", "VALUE")

    clean_price = bond.clean_price_from_ytm(settlement_date, ytm)
    asw = bond.asset_swap_spread(settlement_date, clean_price, flat_curve)
    testCases.print("Discounted on Bond Curve ASW:", asw * 10000)

    # When the libor curve is the Libor curve then the ASW is positive
    libor_curve = build_Ibor_Curve(settlement_date)
    asw = bond.asset_swap_spread(settlement_date, clean_price, libor_curve)
    oas = bond.option_adjusted_spread(
        settlement_date, clean_price, libor_curve)
    testCases.print("Discounted on LIBOR Curve ASW:", asw * 10000)
    testCases.print("Discounted on LIBOR Curve OAS:", oas * 10000)

    p = 90.0
    asw = bond.asset_swap_spread(settlement_date, p, libor_curve)
    oas = bond.option_adjusted_spread(settlement_date, p, libor_curve)
    testCases.print("Deep discount bond at 90 ASW:", asw * 10000)
    testCases.print("Deep discount bond at 90 OAS:", oas * 10000)

    p = 100.0
    asw = bond.asset_swap_spread(settlement_date, p, libor_curve)
    oas = bond.option_adjusted_spread(settlement_date, p, libor_curve)
    testCases.print("Par bond at 100 ASW:", asw * 10000)
    testCases.print("Par bond at 100 OAS:", oas * 10000)

    p = 120.0
    asw = bond.asset_swap_spread(settlement_date, p, libor_curve)
    oas = bond.option_adjusted_spread(settlement_date, p, libor_curve)
    testCases.print("Above par bond at 120 ASW:", asw * 10000)
    testCases.print("Above par bond at 120 OAS:", oas * 10000)

    ##########################################################################
    # https://data.bloomberglp.com/bat/sites/3/2017/07/SF-2017_Paul-Fjeldsted.pdf
    # Page 10 TREASURY NOTE SCREENSHOT
    ##########################################################################

    testCases.banner("BLOOMBERG US TREASURY EXAMPLE")
    settlement_date = Date(21, 7, 2017)
    issue_date = Date(15, 5, 2010)
    maturity_date = Date(15, 5, 2027)
    coupon = 0.02375
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    face = 100.0

    bond = Bond(issue_date,
                maturity_date,
                coupon,
                freq_type,
                accrual_type,
                face)

    testCases.header("FIELD", "VALUE")
    clean_price = 99.7808417

    yld = bond.current_yield(clean_price)
    testCases.print("Current Yield = ", yld)

    ytm = bond.yield_to_maturity(settlement_date, clean_price,
                                 YTMCalcType.UK_DMO)
    testCases.print("UK DMO Yield To Maturity = ", ytm)

    ytm = bond.yield_to_maturity(settlement_date, clean_price,
                                 YTMCalcType.US_STREET)
    testCases.print("US STREET Yield To Maturity = ", ytm)

    ytm = bond.yield_to_maturity(settlement_date, clean_price,
                                 YTMCalcType.US_TREASURY)
    testCases.print("US TREASURY Yield To Maturity = ", ytm)

    full_price = bond.full_price_from_ytm(settlement_date, ytm)
    testCases.print("Full Price = ", full_price)

    clean_price = bond.clean_price_from_ytm(settlement_date, ytm)
    testCases.print("Clean Price = ", clean_price)

    accrued_interest = bond._accrued_interest
    testCases.print("Accrued = ", accrued_interest)

    accddays = bond._accrued_days
    testCases.print("Accrued Days = ", accddays)

    duration = bond.dollar_duration(settlement_date, ytm)
    testCases.print("Dollar Duration = ", duration)

    modified_duration = bond.modified_duration(settlement_date, ytm)
    testCases.print("Modified Duration = ", modified_duration)

    macauley_duration = bond.macauley_duration(settlement_date, ytm)
    testCases.print("Macauley Duration = ", macauley_duration)

    conv = bond.convexity_from_ytm(settlement_date, ytm)
    testCases.print("Convexity = ", conv)

    ##########################################################################
    # Page 11 APPLE NOTE SCREENSHOT
    ##########################################################################

    testCases.banner("BLOOMBERG APPLE CORP BOND EXAMPLE")
    settlement_date = Date(21, 7, 2017)
    issue_date = Date(13, 5, 2012)
    maturity_date = Date(13, 5, 2022)
    coupon = 0.027
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.THIRTY_E_360_ISDA
    face = 100.0

    bond = Bond(issue_date, maturity_date,
                coupon, freq_type, accrual_type, face)

    testCases.header("FIELD", "VALUE")
    clean_price = 101.581564

    yld = bond.current_yield(clean_price)
    testCases.print("Current Yield", yld)

    ytm = bond.yield_to_maturity(settlement_date, clean_price,
                                 YTMCalcType.UK_DMO)
    testCases.print("UK DMO Yield To Maturity", ytm)

    ytm = bond.yield_to_maturity(settlement_date, clean_price,
                                 YTMCalcType.US_STREET)
    testCases.print("US STREET Yield To Maturity", ytm)

    ytm = bond.yield_to_maturity(settlement_date, clean_price,
                                 YTMCalcType.US_TREASURY)
    testCases.print("US TREASURY Yield To Maturity", ytm)

    full_price = bond.full_price_from_ytm(settlement_date, ytm)
    testCases.print("Full Price", full_price)

    clean_price = bond.clean_price_from_ytm(settlement_date, ytm)
    testCases.print("Clean Price", clean_price)

    accddays = bond._accrued_days
    testCases.print("Accrued Days", accddays)

    accrued_interest = bond._accrued_interest
    testCases.print("Accrued", accrued_interest)

    duration = bond.dollar_duration(settlement_date, ytm)
    testCases.print("Dollar Duration", duration)

    modified_duration = bond.modified_duration(settlement_date, ytm)
    testCases.print("Modified Duration", modified_duration)

    macauley_duration = bond.macauley_duration(settlement_date, ytm)
    testCases.print("Macauley Duration", macauley_duration)

    conv = bond.convexity_from_ytm(settlement_date, ytm)
    testCases.print("Convexity", conv)


###############################################################################


def test_BondExDividend():
    issue_date = Date(7, 9, 2000)
    maturity_date = Date(7, 9, 2020)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    face = 100.0
    exDivDays = 7
    testCases.header("LABEL", "VALUE")

    calendar_type = CalendarTypes.UNITED_KINGDOM
    bond = Bond(issue_date, maturity_date, coupon,
                freq_type, accrual_type, face)
    settlement_date = Date(7, 9, 2003)
    accrued = bond.calc_accrued_interest(
        settlement_date, exDivDays, calendar_type)
    testCases.print("SettlementDate:", settlement_date)
    testCases.print("Accrued:", accrued)

    ###########################################################################
    testCases.banner("=======================================================")
    testCases.header("SETTLEMENT", "ACCRUED")

    issue_date = Date(7, 9, 2000)
    maturity_date = Date(7, 9, 2020)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    face = 100.0
    exDivDays = 7

    calendar_type = CalendarTypes.UNITED_KINGDOM
    bond = Bond(issue_date, maturity_date, coupon,
                freq_type, accrual_type, face)

    settlement_date = Date(25, 8, 2010)

    for _ in range(0, 13):
        settlement_date = settlement_date.add_days(1)
        accrued = bond.calc_accrued_interest(
            settlement_date, exDivDays, calendar_type)
        testCases.print(settlement_date, accrued)


###############################################################################


test_Bond()
test_BondExDividend()
testCases.compareTestCases()
