##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import os
import datetime as dt
import pandas as pd
import numpy as np

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


def test_bond():

    import pandas as pd
    path = os.path.join(os.path.dirname(__file__), './data/giltBondPrices.txt')
    bondDataFrame = pd.read_csv(path, sep='\t')
    bondDataFrame['mid'] = 0.5 * (bondDataFrame['bid'] + bondDataFrame['ask'])

    freq_type = FrequencyTypes.SEMI_ANNUAL
    settlement_date = Date(19, 9, 2012)
    face = ONE_MILLION
    ex_div_days = 0

    for accrual_type in DayCountTypes:
        if accrual_type == DayCountTypes.ZERO:
            continue
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
                        coupon, freq_type, accrual_type, ex_div_days)

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
    ex_div_days = 0
    face = 1000000.0

    freq_type = FrequencyTypes.SEMI_ANNUAL

    bond = Bond(issue_date, maturity_date,
                coupon, freq_type, accrualConvention, ex_div_days)

    testCases.header("FIELD", "VALUE")
    dirty_price = bond.dirty_price_from_ytm(settlement_date, y)
    testCases.print("Dirty Price = ", dirty_price)

    clean_price = bond.clean_price_from_ytm(settlement_date, y)
    testCases.print("Clean Price = ", clean_price)

    accrued_interest = bond.accrued_interest(settlement_date, face)
    testCases.print("Accrued = ", accrued_interest)

    ytm = bond.yield_to_maturity(settlement_date, clean_price)
    testCases.print("Yield to Maturity = ", ytm)

    bump = 1e-4
    priceBumpedUp = bond.dirty_price_from_ytm(settlement_date, y + bump)
    testCases.print("Price Bumped Up:", priceBumpedUp)

    priceBumpedDn = bond.dirty_price_from_ytm(settlement_date, y - bump)
    testCases.print("Price Bumped Dn:", priceBumpedDn)

    durationByBump = -(priceBumpedUp - dirty_price) / bump
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
    ex_div_days = 0
    face = 1000000.0

    bond = Bond(issue_date,
                maturity_date,
                coupon,
                freq_type,
                accrual_type,
                ex_div_days,
                CalendarTypes.UNITED_STATES)

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

    dirty_price = bond.dirty_price_from_ytm(settlement_date, ytm,
                                          YTMCalcType.US_TREASURY)

    testCases.print("Dirty Price = ", dirty_price)

    clean_price = bond.clean_price_from_ytm(settlement_date, ytm,
                                            YTMCalcType.US_TREASURY)
    testCases.print("Clean Price = ", clean_price)

    accrued_interest = bond.accrued_interest(settlement_date, face)
    testCases.print("Accrued = ", accrued_interest)

    accddays = bond._accrued_days
    testCases.print("Accrued Days = ", accddays)

    duration = bond.dollar_duration(settlement_date, ytm, YTMCalcType.US_STREET)
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
    ex_div_days = 0
    face = 100.0

    bond = Bond(issue_date, maturity_date,
                coupon, freq_type, accrual_type, ex_div_days)

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

    dirty_price = bond.dirty_price_from_ytm(settlement_date, ytm)
    testCases.print("Dirty Price", dirty_price)

    clean_price = bond.clean_price_from_ytm(settlement_date, ytm)
    testCases.print("Clean Price", clean_price)

    accddays = bond._accrued_days
    testCases.print("Accrued Days", accddays)

    accrued_interest = bond.accrued_interest(settlement_date, face)
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
    ex_div_days = 7
    testCases.header("LABEL", "VALUE")

    bond = Bond(issue_date, maturity_date, coupon,
                freq_type, accrual_type, ex_div_days)
    settlement_date = Date(7, 9, 2003)
    accrued = bond.accrued_interest(settlement_date, face)

    testCases.print("SettlementDate:", settlement_date)
    testCases.print("Accrued:", accrued)

    ###########################################################################
    testCases.banner("=======================================================")
    testCases.header("SETTLEMENT", "DIRTY PRICE", "ACCRUED", "CLEAN PRICE")

    issue_date = Date(7, 9, 2000)
    maturity_date = Date(7, 9, 2020)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    face = 100.0
    ex_div_days = 7

    bond = Bond(issue_date, maturity_date, coupon,
                freq_type, accrual_type, ex_div_days)

    settlement_date = Date(25, 8, 2010)

    ytm = 0.05

    for _ in range(0, 13):
        settlement_date = settlement_date.add_days(1)
        accrued = bond.accrued_interest(
            settlement_date, face)
        dirty_price = bond.dirty_price_from_ytm(
            settlement_date, ytm)
        clean_price = dirty_price - accrued
        testCases.print(settlement_date, dirty_price, accrued, clean_price)

###############################################################################


def test_BondPaymentDates():

    from financepy.products.bonds.bond import Bond
    from financepy.utils import Date, DayCountTypes, FrequencyTypes

    bond = Bond(
            issue_date=Date(7, 6, 2021),
            maturity_date=Date(7, 6, 2031),
            coupon=0.0341,
            freq_type=FrequencyTypes.ANNUAL,
            accrual_type=DayCountTypes.ACT_ACT_ISDA
    )
    bond._calculate_payment_dates()

    if 1 == 0:
        print(bond._flow_amounts)
        print(bond._coupon_dates)
        print(bond._payment_dates)

###############################################################################


def test_Bond_ror():

    test_case_file = 'test_cases_bond_ror.csv'
    df = pd.read_csv('./data/' + test_case_file, parse_dates=['buy_date', 'sell_date'])
    # A 10-year bond with 1 coupon per year. code: 210215
    bond = Bond(
        issue_date=Date(13, 9, 2021),
        maturity_date=Date(13, 9, 2031),
        coupon=0.0312,
        freq_type=FrequencyTypes.ANNUAL,
        accrual_type=DayCountTypes.ACT_ACT_ICMA
    )
    testCases.header('bond_code', 'buy_date', 'buy_ytm', 'buy_price', 'sell_date', 'sell_ytm', 'sell_price',
                     'simple_return', 'irr')
    for row in df.itertuples(index=False):
        buy_date = Date(row.buy_date.day, row.buy_date.month, row.buy_date.year)
        sell_date = Date(row.sell_date.day, row.sell_date.month, row.sell_date.year)
        buy_price = bond.dirty_price_from_ytm(buy_date, row.buy_ytm, YTMCalcType.US_STREET)
        sell_price = bond.dirty_price_from_ytm(sell_date, row.sell_ytm, YTMCalcType.US_STREET)
        simple, irr, pnl = bond.calc_ror(buy_date, sell_date, row.buy_ytm, row.sell_ytm)
        testCases.print(row.bond_code, buy_date, row.buy_ytm, buy_price, sell_date, row.sell_ytm, sell_price,
                        simple, irr)


###############################################################################

def test_Bond_eom():

    # Bonds that mature on an EOM date have flows on EOM dates
    issue_date = Date(30, 11, 2022)
    settle_date = Date(6, 2, 2023)
    maturity_date = Date(30, 11, 2024)
    coupon = 0.045
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    ex_div_days = 0

    bond = Bond(issue_date, maturity_date, coupon, freq_type, accrual_type, ex_div_days)

    ai = bond.accrued_interest(settle_date) # should be 8406.593406

###############################################################################

def test_key_rate_durations():

#    print("Complete test case logging")

    issue_date = Date(31, 7, 2022)
    maturity_date = Date(31, 7, 2027)
    coupon = 0.0275
    ex_div_days = 0

    accrual_type, freq_type, settlementDays, exDiv, calendar = get_bond_market_conventions(
    BondMarkets.UNITED_STATES)

    bond = Bond(issue_date, maturity_date, coupon,
                freq_type, accrual_type, ex_div_days)

    settlement_date = Date(24, 4, 2023)

    ytm = 3.725060/100

    key_rate_tenors, key_rate_durations = bond.key_rate_durations(
        settlement_date, ytm)

#    print(key_rate_tenors)
#    print(key_rate_durations)

###############################################################################


def test_key_rate_durations_Bloomberg_example():

    accrual_type, frequencyType, settlementDays, exDiv, calendar = \
    get_bond_market_conventions(BondMarkets.UNITED_STATES)

    # interest accrues on this date. Issue date is 01/08/2022
    issue_date = Date(31, 7, 2022)
    maturity_date = Date(31, 7, 2027)
    coupon = 2.75/100.0
    ex_div_days = 0

    accrual_type, freq_type, settlementDays, exDiv, calendar = get_bond_market_conventions(
    BondMarkets.UNITED_STATES)

    bond = Bond(issue_date, maturity_date, coupon,
                freq_type, accrual_type, ex_div_days)

    settlement_date = Date(24, 4, 2023)

    # US Street yield on Bloomberg as of 20 April 2023
    # with settle date 24 April 2023
    ytm = 3.725060/100

    # Details of yields of market bonds at KRD maturity points
    my_tenors = np.array([0.5,  1,  2,  3,  5,  7,  10])
    my_rates = np.array([5.0367, 4.7327, 4.1445, 3.8575, 3.6272,  3.5825,  3.5347])/100

    key_rate_tenors, key_rate_durations = bond.key_rate_durations(settlement_date,
                                                                  ytm,
                                                                  key_rate_tenors = my_tenors,
                                                                  rates = my_rates)

#    print(key_rate_tenors)
#    print(key_rate_durations)

    # Differences due to bonds not sitting exactly on these maturity points ? Did BBG interpolate ?

###############################################################################

from financepy.market.curves.discount_curve_flat import DiscountCurveFlat

def test_oas():

    issue_date = Date(15, 5, 2010)
    maturity_date = Date(15, 5, 2027)
    coupon = 0.02375
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    face = 100.0 # By setting the face to 100 we expect a price of par to be 100.0

    bond = Bond(issue_date, maturity_date, coupon, freq_type, accrual_type)

    liborFlatRate = 0.0275
    settlement_date = Date(21, 7, 2017)

    liborFlatCurve = DiscountCurveFlat(settlement_date, liborFlatRate, FrequencyTypes.SEMI_ANNUAL)

    clean_price = 99.780842  # I specified face to be 100 - if face is 1 then this must be 0.99780842

    oas = bond.option_adjusted_spread(settlement_date, clean_price, liborFlatCurve) * 10000

    if (oas - (-34.95)) > 0.01:
        print("OAS incorrect")

###############################################################################

def test_div_dates():

    issueDate = Date(15, 5, 2020)
    maturityDate = Date(15, 5, 2035)
    coupon = 0.02375
    freqType = FrequencyTypes.SEMI_ANNUAL
    accrualType = DayCountTypes.ACT_ACT_ICMA
    face = 125000
    ex_div_days = 10

    bond = Bond(issueDate, maturityDate, coupon, freqType, accrualType, ex_div_days)

    print(bond)

    cleanPrice = 99.7808417  # if face is 1 then this must be 0.99780842

    settlementDate = Date(15, 5, 2023)
    print(bond.bond_payments(settlementDate, face))

    current_yield = bond.current_yield(cleanPrice)*100
    print("Currnt Yield: %10.5f %%" % (current_yield))

    ytm = bond.yield_to_maturity(settlementDate, cleanPrice) * 100.0
    print("Yield to Mat: %10.5f %%" % (ytm))

###############################################################################

test_bond()
test_oas()
test_BondExDividend()
test_BondPaymentDates()
test_Bond_ror()
test_Bond_eom()
test_key_rate_durations()
test_key_rate_durations_Bloomberg_example()

testCases.compareTestCases()
