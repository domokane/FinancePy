###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import matplotlib.pyplot as plt
import numpy as np
import time as time

import sys
sys.path.append("..")

from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.calendar import CalendarTypes
from financepy.products.rates.FinIborFRA import FinIborFRA
from financepy.products.rates.FinIborFuture import FinIborFuture
from financepy.products.rates.FinOIS import FinOIS
from financepy.products.rates.FinOISCurve import OISCurve
from financepy.products.rates.FinIborDeposit import FinIborDeposit
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.market.curves.interpolator import FinInterpTypes
from financepy.utils.FinGlobalTypes import FinSwapTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

PLOT_GRAPHS = False

###############################################################################


def test_FinOISFRAsOnly():

    # TO DO FIX THIS
    valuation_date = Date(23, 2, 2018)

    spotDays = 0
    settleDt = valuation_date.addWeekDays(spotDays)

    depoDCCType = DayCountTypes.ACT_360
    notional = 100.0

    payFixed = True

    calendar_type = CalendarTypes.TARGET
    fras = []

    # 1 x 4 FRA
    fraRate = 0.04
    frasettleDt = settleDt.addMonths(1)
    fraMaturityDate = settleDt.addMonths(4)
    fra = FinIborFRA(frasettleDt, fraMaturityDate, fraRate,
                      depoDCCType, notional, payFixed, calendar_type)
    fras.append(fra)

    # 4 x 7 FRA
    fraRate = 0.08
    frasettleDt = settleDt.addMonths(4)
    fraMaturityDate = settleDt.addMonths(7)
    fra = FinIborFRA(frasettleDt, fraMaturityDate, fraRate,
                      depoDCCType, notional, payFixed, calendar_type)
    fras.append(fra)

    swaps = []

    libor_curve = OISCurve(settleDt,
                           fras,
                           swaps)

    testCases.header("DATE", "MATDATE", "VALUE")

    """ Check calibration """
    for fra in fras:
        v = fra.value(settleDt, libor_curve)
        testCases.print("FRA:", fra._maturity_date, v)

###############################################################################


def test_FinOISDepositsFRAsSwaps():

    valuation_date = Date(18, 9, 2019)

    dccType = DayCountTypes.THIRTY_E_360_ISDA
    depos = []

    spotDays = 0
    settleDt = valuation_date.addWeekDays(spotDays)

    depoDCCType = DayCountTypes.ACT_360
    notional = 100.0
    calendar_type = CalendarTypes.TARGET
    depos = []

    # 1 month
    deposit_rate = 0.04
    maturity_date = settleDt.addMonths(1)
    depo = FinIborDeposit(settleDt, maturity_date, deposit_rate,
                          depoDCCType, notional, calendar_type)
    depos.append(depo)
    
    fras = []
    # 1 x 4 FRA
    fraRate = 0.04
    frasettleDt = settleDt.addMonths(9)
    fraMaturityDate = settleDt.addMonths(13)
    fra = FinIborFRA(frasettleDt, fraMaturityDate, fraRate, dccType)
    fras.append(fra)

    # 4 x 7 FRA
    fraRate = 0.03
    frasettleDt = settleDt.addMonths(13)
    fraMaturityDate = settleDt.addMonths(17)
    fra = FinIborFRA(frasettleDt, fraMaturityDate, fraRate, dccType)
    fras.append(fra)

    # 4 x 7 FRA
    fraRate = 0.07
    frasettleDt = settleDt.addMonths(17)
    fraMaturityDate = settleDt.addMonths(21)
    fra = FinIborFRA(frasettleDt, fraMaturityDate, fraRate, dccType)
    fras.append(fra)

    swaps = []
    fixedDCCType = DayCountTypes.ACT_365F
    fixedFreqType = FrequencyTypes.SEMI_ANNUAL

    swap_rate = 0.05
#    maturity_date = settleDt.addMonths(24)
#    swap = FinIborSwap(settleDt, maturity_date, swap_rate, fixedFreqType,
#                        fixedDCCType)
#    swaps.append(swap)

    fixed_legType = Finfixed_legTypes.PAY
    maturity_date = settleDt.addMonths(36)
    swap = FinOIS(settleDt, maturity_date, fixed_legType, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settleDt.addMonths(48)
    swap = FinOIS(settleDt, maturity_date, fixed_legType, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settleDt.addMonths(60)
    swap = FinOIS(settleDt, maturity_date, fixed_legType, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settleDt.addMonths(72)
    swap = FinOIS(settleDt, maturity_date, fixed_legType, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settleDt.addMonths(84)
    swap = FinOIS(settleDt, maturity_date, fixed_legType, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settleDt.addMonths(96)
    swap = FinOIS(settleDt, maturity_date, fixed_legType, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settleDt.addMonths(108)
    swap = FinOIS(settleDt, maturity_date, fixed_legType, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settleDt.addMonths(120)
    swap = FinOIS(settleDt, maturity_date, fixed_legType, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settleDt.addMonths(132)
    swap = FinOIS(settleDt, maturity_date, fixed_legType, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settleDt.addMonths(144)
    swap = FinOIS(settleDt, maturity_date, fixed_legType, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settleDt.addMonths(180)
    swap = FinOIS(settleDt, maturity_date, fixed_legType, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settleDt.addMonths(240)
    swap = FinOIS(settleDt, maturity_date, fixed_legType, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settleDt.addMonths(300)
    swap = FinOIS(settleDt, maturity_date, fixed_legType, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settleDt.addMonths(360)
    swap = FinOIS(settleDt, maturity_date, fixed_legType, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    libor_curve = OISCurve(valuation_date,
                           depos,
                           fras,
                           swaps)

    df = libor_curve.df(settleDt)

    testCases.header("SETTLEMENT DATE", "DF")
    testCases.print(str(settleDt), df)
    testCases.header("DATE", "DF")

    for deposit in depos:
        df = libor_curve.df(deposit._maturity_date)
        testCases.print(str(deposit._maturity_date), df)

    for swap in swaps:
        df = libor_curve.df(swap._maturity_date)
        testCases.print(str(swap._maturity_date), df)


###############################################################################
###############################################################################
# https://github.com/vilen22/curve-building/blob/master/Bloomberg%20Curve%20Building%20Replication.xlsx
###############################################################################
###############################################################################
# AGREEEMENT IS VERY CLOSE - NOT SURE ABOUT SIZE OF LAST PAYMENT ON FIXED LEG!


def futureToFRARate(price, convexity):
    futRate = (100-price)/100
    if convexity < 0:
        fraRate = futRate + convexity/100.0
    else:
        fraRate = futRate - convexity/100.0

    return fraRate

###############################################################################


def test_FinOISDepositsFuturesSwaps():

    
    spotDate = Date(6, 6, 2018)
    spotDays = 0
    settleDt = spotDate.addWeekDays(spotDays)
    depoDCCType = DayCountTypes.THIRTY_E_360_ISDA

    depo = FinIborDeposit(settleDt, "1D", 1.712/100.0, depoDCCType)
    depos = [depo]

    fras = []

    fraRate = futureToFRARate(97.6675, -0.00005)
    frasettleDt = spotDate.nextIMMDate()
    fraMaturityDate = frasettleDt.nextIMMDate()
    fra = FinIborFRA(frasettleDt, fraMaturityDate, fraRate, depoDCCType)
    fras.append(fra)

    fraRate = futureToFRARate(97.5200, -0.00060)
    frasettleDt = fraMaturityDate
    fraMaturityDate = frasettleDt.nextIMMDate()
    fra = FinIborFRA(frasettleDt, fraMaturityDate, fraRate, depoDCCType)
    fras.append(fra)

    fraRate = futureToFRARate(97.3550, -0.00146)
    frasettleDt = fraMaturityDate
    fraMaturityDate = frasettleDt.nextIMMDate()
    fra = FinIborFRA(frasettleDt, fraMaturityDate, fraRate, depoDCCType)
    fras.append(fra)

    fraRate = futureToFRARate(97.2450, -0.00263)
    frasettleDt = fraMaturityDate
    fraMaturityDate = frasettleDt.nextIMMDate()
    fra = FinIborFRA(frasettleDt, fraMaturityDate, fraRate, depoDCCType)
    fras.append(fra)

    fraRate = futureToFRARate(97.1450, -0.00411)
    frasettleDt = fraMaturityDate
    fraMaturityDate = frasettleDt.nextIMMDate()
    fra = FinIborFRA(frasettleDt, fraMaturityDate, fraRate, depoDCCType)
    fras.append(fra)

    fraRate = futureToFRARate(97.0750, -0.00589)
    frasettleDt = frasettleDt.nextIMMDate()
    fraMaturityDate = frasettleDt.nextIMMDate()
    fra = FinIborFRA(frasettleDt, fraMaturityDate, fraRate, depoDCCType)
    fras.append(fra)

    ###########################################################################

    spotDays = 2
    start_date = spotDate.addWeekDays(spotDays)

    swaps = []
    fixed_legType = FinSwapTypes.PAY
    fixedDCCType = DayCountTypes.THIRTY_E_360
    fixedFreqType = FrequencyTypes.SEMI_ANNUAL
    floatFreqType = FrequencyTypes.QUARTERLY
    notional = 1000000
    floatSpread = 0.0
    floatDCCType = DayCountTypes.ACT_360
    calendar_type = CalendarTypes.US
    busDayAdjustRule = BusDayAdjustTypes.PRECEDING

    swap_rate = 0.02776305
    payment_lag = 1

    swap = FinOIS(start_date, "2Y", fixed_legType,
                  swap_rate, fixedFreqType, fixedDCCType, notional,
                  payment_lag, floatSpread, floatFreqType, floatDCCType,
                  calendar_type, busDayAdjustRule)

    swaps.append(swap)

    libor_curve = OISCurve(spotDate, depos, fras, swaps)

    times = np.linspace(0.0, 2.0, 25)
    dates = spotDate.addYears(times)
    zeroRates = libor_curve.zeroRate(dates)
    fwd_rates = libor_curve.fwd(dates)

    if PLOT_GRAPHS:
        plt.figure(figsize=(8, 6))
        plt.plot(times, zeroRates*100, label="zero rates")
        plt.plot(times, fwd_rates*100, label="fwd rates")
        plt.xlabel("Times")
        plt.ylabel("CC forward rates")
        plt.legend()

        print("==============================================================")
        for fra in fras:
            print(fra)
        print("==============================================================")

        end_date = spotDate
        df = libor_curve.df(end_date)
        print(end_date, df)

        end_date = settleDt
        df = libor_curve.df(end_date)
        print(end_date, df)

        end_date = Date(20, 6, 2018)
        df = libor_curve.df(end_date)
        print(end_date, df)

        for fra in fras:
            end_date = fra._maturity_date
            df = libor_curve.df(end_date)
            print(end_date, df)

        for swap in swaps:
            end_date = swap._maturity_date
            df = libor_curve.df(end_date)
            print(end_date, df)

        swap.printFixedLegPV(spotDate)
        swap.printFloatLegPV(spotDate)

#        print(libor_curve)

###############################################################################


def test_derivativePricingExample():

    valuation_date = Date(10, 11, 2011)

    # We do the O/N rate which settles on trade date
    spotDays = 0
    settleDt = valuation_date.addWeekDays(spotDays)

    fras = []

    swaps = []
    day_count_type = DayCountTypes.THIRTY_E_360_ISDA
#    day_count_type = DayCountTypes.ACT_360
    freq_type = FrequencyTypes.SEMI_ANNUAL
    fixed_legType = Finfixed_legTypes.PAY
    
    swap_rate = 0.0058
    swap = FinOIS(settleDt, "1Y", fixed_legType, swap_rate, freq_type, day_count_type)
    swaps.append(swap)

    swap_rate = 0.0060
    swap = FinOIS(settleDt, "2Y", fixed_legType, swap_rate, freq_type, day_count_type)
    swaps.append(swap)

    swap_rate = 0.0072
    swap = FinOIS(settleDt, "3Y", fixed_legType, swap_rate, freq_type, day_count_type)
    swaps.append(swap)

    swap_rate = 0.0096
    swap = FinOIS(settleDt, "4Y", fixed_legType, swap_rate, freq_type, day_count_type)
    swaps.append(swap)

    swap_rate = 0.0124
    swap = FinOIS(settleDt, "5Y", fixed_legType, swap_rate, freq_type, day_count_type)
    swaps.append(swap)

    swap_rate = 0.0173
    swap = FinOIS(settleDt, "7Y", fixed_legType, swap_rate, freq_type, day_count_type)
    swaps.append(swap)

    swap_rate = 0.0219
    swap = FinOIS(settleDt, "10Y", fixed_legType, swap_rate, freq_type, day_count_type)
    swaps.append(swap)

    swap_rate = 0.0283
    swap = FinOIS(settleDt, "30Y", fixed_legType, swap_rate, freq_type, day_count_type)
    swaps.append(swap)

    numRepeats = 10
    start = time.time()

    for _ in range(0, numRepeats):
        _ = OISCurve(valuation_date, fras, swaps,
                     FinInterpTypes.FLAT_FWD_RATES)

    end = time.time()
    elapsed1 = end - start

    start = time.time()

    for _ in range(0, numRepeats):
        _ = OISCurve(valuation_date, fras, swaps,
                     FinInterpTypes.LINEAR_SWAP_RATES)

    end = time.time()
    elapsed2 = end - start

    testCases.header("METHOD", "TIME")
    testCases.print("NON-LINEAR SOLVER BOOTSTRAP", elapsed1/numRepeats)
    testCases.print("LINEAR SWAP BOOTSTRAP", elapsed2/numRepeats)

###############################################################################


def test_bloombergPricingExample():

    """ This is an example of a replication of a BBG example from
    https://github.com/vilen22/curve-building/blob/master/Bloomberg%20Curve%20Building%20Replication.xlsx
    """

    valuation_date = Date(6, 6, 2018)

    # We do the O/N rate which settles on trade date
    spotDays = 0
    settleDt = valuation_date.addWeekDays(spotDays)
    accrual = DayCountTypes.THIRTY_E_360

    depo = FinIborDeposit(settleDt, "1D", 1.712/100.0, accrual)
    depos = [depo]
    
    futs = []
    fut = FinIborFuture(valuation_date, 1); futs.append(fut)
    fut = FinIborFuture(valuation_date, 2); futs.append(fut)
    fut = FinIborFuture(valuation_date, 3); futs.append(fut)
    fut = FinIborFuture(valuation_date, 4); futs.append(fut)
    fut = FinIborFuture(valuation_date, 5); futs.append(fut)
    fut = FinIborFuture(valuation_date, 6); futs.append(fut)

    fras = [None]*6
    fras[0] = futs[0].toFRA(97.6675, -0.00005)
    fras[1] = futs[1].toFRA(97.5200, -0.00060)
    fras[2] = futs[2].toFRA(97.3550, -0.00146)
    fras[3] = futs[3].toFRA(97.2450, -0.00263)
    fras[4] = futs[4].toFRA(97.1450, -0.00411)
    fras[5] = futs[5].toFRA(97.0750, -0.00589)

    accrual = DayCountTypes.THIRTY_E_360
    freq = FrequencyTypes.SEMI_ANNUAL
    spotDays = 2
    settleDt = valuation_date.addWeekDays(spotDays)
    payRec = FinSwapTypes.PAY
    lag = 1 # Not used

    swaps = []
    swap = FinOIS(settleDt, "2Y", payRec, (2.77417+2.77844)/200, freq, accrual); swaps.append(swap)
    swap = FinOIS(settleDt, "3Y", payRec, (2.86098+2.86582)/200, freq, accrual); swaps.append(swap)
    swap = FinOIS(settleDt, "4Y", payRec, (2.90240+2.90620)/200, freq, accrual); swaps.append(swap)
    swap = FinOIS(settleDt, "5Y", payRec, (2.92944+2.92906)/200, freq, accrual); swaps.append(swap)
    swap = FinOIS(settleDt, "6Y", payRec, (2.94001+2.94499)/200, freq, accrual); swaps.append(swap)
    swap = FinOIS(settleDt, "7Y", payRec, (2.95352+2.95998)/200, freq, accrual); swaps.append(swap)
    swap = FinOIS(settleDt, "8Y", payRec, (2.96830+2.97400)/200, freq, accrual); swaps.append(swap)
    swap = FinOIS(settleDt, "9Y", payRec, (2.98403+2.98817)/200, freq, accrual); swaps.append(swap)
    swap = FinOIS(settleDt, "10Y", payRec, (2.99716+3.00394)/200, freq, accrual); swaps.append(swap)
    swap = FinOIS(settleDt, "11Y", payRec, (3.01344+3.01596)/200, freq, accrual); swaps.append(swap)
    swap = FinOIS(settleDt, "12Y", payRec, (3.02276+3.02684)/200, freq, accrual); swaps.append(swap)
    swap = FinOIS(settleDt, "15Y", payRec, (3.04092+3.04508)/200, freq, accrual); swaps.append(swap)
    swap = FinOIS(settleDt, "20Y", payRec, (3.04417+3.05183)/200, freq, accrual); swaps.append(swap)
    swap = FinOIS(settleDt, "25Y", payRec, (3.03219+3.03621)/200, freq, accrual); swaps.append(swap)
    swap = FinOIS(settleDt, "30Y", payRec, (3.01030+3.01370)/200, freq, accrual); swaps.append(swap)
    swap = FinOIS(settleDt, "40Y", payRec, (2.96946+2.97354)/200, freq, accrual); swaps.append(swap)
    swap = FinOIS(settleDt, "50Y", payRec, (2.91552+2.93748)/200, freq, accrual); swaps.append(swap)

    oisCurve = OISCurve(valuation_date, depos, fras, swaps)

#    swaps[0]._fixed_leg.printValuation()
#    swaps[0]._floatLeg.printValuation()
    
    # The valuation of 53714.55 is very close to the spreadsheet value 53713.96
    principal = 0.0

    testCases.header("VALUATION TO TODAY DATE"," PV")
    testCases.print("VALUE:", swaps[0].value(valuation_date, oisCurve, None))
    testCases.print("FIXED:", -swaps[0]._fixed_leg.value(valuation_date, oisCurve))
    testCases.print("FLOAT:", swaps[0]._floatLeg.value(valuation_date, oisCurve, None))

    testCases.header("VALUATION TO SWAP SETTLEMENT DATE"," PV")
    testCases.print("VALUE:", swaps[0].value(settleDt, oisCurve, None))
    testCases.print("FIXED:", -swaps[0]._fixed_leg.value(settleDt, oisCurve))
    testCases.print("FLOAT:", swaps[0]._floatLeg.value(settleDt, oisCurve, None))

    # swaps[0].printFixedLegPV()
    # swaps[0].printFloatLegPV()

###############################################################################


test_bloombergPricingExample()
#test_derivativePricingExample()
#test_FinOISFRAsOnly()
#test_FinOISDepositsFRAsSwaps()
#test_FinOISDepositsFuturesSwaps()

testCases.compareTestCases()
