###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
import matplotlib.pyplot as plt
import numpy as np
import time as time

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinCalendar import FinCalendarTypes
from financepy.products.libor.FinLiborCurve import FinLiborCurve
from financepy.products.libor.FinLiborFRA import FinLiborFRA
from financepy.products.libor.FinLiborFuture import FinLiborFuture
from financepy.products.libor.FinLiborDeposit import FinLiborDeposit
from financepy.products.libor.FinLiborSwap import FinLiborSwap
from financepy.finutils.FinCalendar import FinBusDayAdjustTypes
from financepy.market.curves.FinInterpolate import FinInterpTypes
from financepy.finutils.FinMath import ONE_MILLION

sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)

PLOT_GRAPHS = False

###############################################################################


def test_FinLiborDepositsOnly():

    # I have used the following useful blog post by Ioannis Rigopoulos for this
    # https://blog.deriscope.com/index.php/en/yield-curve-excel-quantlib-deposit

    valuationDate = FinDate(2018, 2, 23)

    spotDays = 0
    settlementDate = valuationDate.addWorkDays(spotDays)

    depoDCCType = FinDayCountTypes.ACT_360
    notional = 100.0
    calendarType = FinCalendarTypes.TARGET
    depos = []

    # 1 month
    depositRate = 0.04
    maturityDate = settlementDate.addMonths(1)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate,
                           depoDCCType, notional, calendarType)
    depos.append(depo)

    # 2 months
    depositRate = 0.04
    maturityDate = settlementDate.addMonths(2)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate,
                           depoDCCType, notional, calendarType)
    depos.append(depo)

    # 6 months
    depositRate = 0.04
    maturityDate = settlementDate.addMonths(6)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate,
                           depoDCCType, notional, calendarType)
    depos.append(depo)

    # 1 year
    depositRate = 0.04
    maturityDate = settlementDate.addMonths(12)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate,
                           depoDCCType, notional, calendarType)
    depos.append(depo)

    fras = []
    swaps = []

    liborCurve = FinLiborCurve(settlementDate,
                               depos,
                               fras,
                               swaps)

    testCases.header("LABEL", "DATE", "VALUE")

    ''' Check calibration '''
    for depo in depos:
        v = depo.value(settlementDate, liborCurve)
        testCases.print("DEPO", depo._maturityDate, v)

###############################################################################


def test_FinLiborFRAsOnly():

    # TO DO FIX THIS
    valuationDate = FinDate(2018, 2, 23)

    spotDays = 0
    settlementDate = valuationDate.addWorkDays(spotDays)

    depoDCCType = FinDayCountTypes.ACT_360
    notional = 100.0

    payFixed = True

    calendarType = FinCalendarTypes.TARGET
    fras = []

    # 1 x 4 FRA
    fraRate = 0.04
    fraSettlementDate = settlementDate.addMonths(1)
    fraMaturityDate = settlementDate.addMonths(4)
    fra = FinLiborFRA(fraSettlementDate, fraMaturityDate, fraRate,
                      depoDCCType, notional, payFixed, calendarType)
    fras.append(fra)

    # 4 x 7 FRA
    fraRate = 0.08
    fraSettlementDate = settlementDate.addMonths(4)
    fraMaturityDate = settlementDate.addMonths(7)
    fra = FinLiborFRA(fraSettlementDate, fraMaturityDate, fraRate,
                      depoDCCType, notional, payFixed, calendarType)
    fras.append(fra)

    depos = []
    swaps = []

    liborCurve = FinLiborCurve(settlementDate,
                               depos,
                               fras,
                               swaps)

    testCases.header("DATE", "MATDATE", "VALUE")

    ''' Check calibration '''
    for fra in fras:
        v = fra.value(settlementDate, liborCurve)
        testCases.print("FRA:", fra._maturityDate, v)

###############################################################################


def test_FinLiborDepositsFRAsSwaps():

    valuationDate = FinDate(2019, 9, 18)

    dccType = FinDayCountTypes.THIRTY_E_360_ISDA
    depos = []

    spotDays = 0
    settlementDate = valuationDate.addWorkDays(spotDays)

    depositRate = 0.050
    maturityDate = settlementDate.addMonths(1)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate, dccType)
    depos.append(depo)

    maturityDate = settlementDate.addMonths(2)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate, dccType)
    depos.append(depo)

    maturityDate = settlementDate.addMonths(3)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate, dccType)
    depos.append(depo)

    maturityDate = settlementDate.addMonths(6)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate, dccType)
    depos.append(depo)

    maturityDate = settlementDate.addMonths(9)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate, dccType)
    depos.append(depo)

    maturityDate = settlementDate.addMonths(12)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate, dccType)
    depos.append(depo)

    fras = []
    # 1 x 4 FRA
    fraRate = 0.04
    fraSettlementDate = settlementDate.addMonths(9)
    fraMaturityDate = settlementDate.addMonths(13)
    fra = FinLiborFRA(fraSettlementDate, fraMaturityDate, fraRate, dccType)
    fras.append(fra)

    # 4 x 7 FRA
    fraRate = 0.03
    fraSettlementDate = settlementDate.addMonths(13)
    fraMaturityDate = settlementDate.addMonths(17)
    fra = FinLiborFRA(fraSettlementDate, fraMaturityDate, fraRate, dccType)
    fras.append(fra)

    # 4 x 7 FRA
    fraRate = 0.07
    fraSettlementDate = settlementDate.addMonths(17)
    fraMaturityDate = settlementDate.addMonths(21)
    fra = FinLiborFRA(fraSettlementDate, fraMaturityDate, fraRate, dccType)
    fras.append(fra)

    swaps = []
    fixedDCCType = FinDayCountTypes.ACT_365F
    fixedFreqType = FinFrequencyTypes.SEMI_ANNUAL

    swapRate = 0.05
#    maturityDate = settlementDate.addMonths(24)
#    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType,
#                        fixedDCCType)
#    swaps.append(swap)

    maturityDate = settlementDate.addMonths(36)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(48)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(60)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(72)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(84)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(96)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(108)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(120)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(132)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(144)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(180)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(240)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(300)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(360)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    liborCurve = FinLiborCurve(valuationDate,
                               depos,
                               fras,
                               swaps)

    df = liborCurve.df(settlementDate)

    testCases.header("SETTLEMENT DATE", "DF")
    testCases.print(str(settlementDate), df)
    testCases.header("DATE", "DF")

    for deposit in depos:
        df = liborCurve.df(deposit._maturityDate)
        testCases.print(str(deposit._maturityDate), df)

    for swap in swaps:
        df = liborCurve.df(swap._maturityDate)
        testCases.print(str(swap._maturityDate), df)


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


def test_FinLiborDepositsFuturesSwaps():

    spotDate = FinDate(6, 6, 2018)
    spotDays = 0
    settlementDate = spotDate.addWorkDays(spotDays)
    depoDCCType = FinDayCountTypes.ACT_360
    depos = []
    depositRate = 0.0231381
    depo = FinLiborDeposit(settlementDate, "3M", depositRate, depoDCCType)
    depos.append(depo)

    depositRate = 0.027
    depo = FinLiborDeposit(settlementDate, "3M", depositRate, depoDCCType)
    depos.append(depo)

    depos = []
    depo = FinLiborDeposit(settlementDate, "1M", 0.0230, depoDCCType)
    depos.append(depo)
    depo = FinLiborDeposit(settlementDate, "2M", 0.0235, depoDCCType)
    depos.append(depo)
    depo = FinLiborDeposit(settlementDate, "3M", 0.0240, depoDCCType)
    depos.append(depo)

    fras = []

    fraRate = futureToFRARate(97.6675, -0.00005)
    fraSettlementDate = spotDate.nextIMMDate()
    fraMaturityDate = fraSettlementDate.nextIMMDate()
    fra = FinLiborFRA(fraSettlementDate, fraMaturityDate, fraRate, depoDCCType)
    fras.append(fra)

    fraRate = futureToFRARate(97.5200, -0.00060)
    fraSettlementDate = fraMaturityDate
    fraMaturityDate = fraSettlementDate.nextIMMDate()
    fra = FinLiborFRA(fraSettlementDate, fraMaturityDate, fraRate, depoDCCType)
    fras.append(fra)

    fraRate = futureToFRARate(97.3550, -0.00146)
    fraSettlementDate = fraMaturityDate
    fraMaturityDate = fraSettlementDate.nextIMMDate()
    fra = FinLiborFRA(fraSettlementDate, fraMaturityDate, fraRate, depoDCCType)
    fras.append(fra)

    fraRate = futureToFRARate(97.2450, -0.00263)
    fraSettlementDate = fraMaturityDate
    fraMaturityDate = fraSettlementDate.nextIMMDate()
    fra = FinLiborFRA(fraSettlementDate, fraMaturityDate, fraRate, depoDCCType)
    fras.append(fra)

    fraRate = futureToFRARate(97.1450, -0.00411)
    fraSettlementDate = fraMaturityDate
    fraMaturityDate = fraSettlementDate.nextIMMDate()
    fra = FinLiborFRA(fraSettlementDate, fraMaturityDate, fraRate, depoDCCType)
    fras.append(fra)

    fraRate = futureToFRARate(97.0750, -0.00589)
    fraSettlementDate = fraSettlementDate.nextIMMDate()
    fraMaturityDate = fraSettlementDate.nextIMMDate()
    fra = FinLiborFRA(fraSettlementDate, fraMaturityDate, fraRate, depoDCCType)
    fras.append(fra)

    ###########################################################################

    spotDays = 2
    startDate = spotDate.addWorkDays(spotDays)

    swaps = []
    fixedDCCType = FinDayCountTypes.THIRTY_E_360
    fixedFreqType = FinFrequencyTypes.SEMI_ANNUAL
    floatFreqType = FinFrequencyTypes.QUARTERLY
    notional = 1000000
    floatSpread = 0.0
    floatDCCType = FinDayCountTypes.ACT_360
    payFixed = True
    calendarType = FinCalendarTypes.US
    busDayAdjustRule = FinBusDayAdjustTypes.PRECEDING

    swapRate = 0.02776305

    swap = FinLiborSwap(startDate, "2Y", swapRate,
                        fixedFreqType, fixedDCCType, notional,
                        floatSpread, floatFreqType, floatDCCType,
                        payFixed, calendarType, busDayAdjustRule)

    swaps.append(swap)

    liborCurve = FinLiborCurve(spotDate, depos, fras, swaps)

    times = np.linspace(0.0, 2.0, 25)
    dates = spotDate.addYears(times)
    zeroRates = liborCurve.zeroRate(dates)
    fwdRates = liborCurve.fwd(dates)

    if PLOT_GRAPHS:
        plt.figure(figsize=(8, 6))
        plt.plot(times, zeroRates*100, label="zero rates")
        plt.plot(times, fwdRates*100, label="fwd rates")
        plt.xlabel("Times")
        plt.ylabel("CC forward rates")
        plt.legend()

        print("==============================================================")
        for fra in fras:
            print(fra)
        print("==============================================================")

        endDate = spotDate
        df = liborCurve.df(endDate)
        print(endDate, df)

        endDate = settlementDate
        df = liborCurve.df(endDate)
        print(endDate, df)

        endDate = FinDate(20, 6, 2018)
        df = liborCurve.df(endDate)
        print(endDate, df)

        for depo in depos:
            endDate = depo._maturityDate
            df = liborCurve.df(endDate)
            print(endDate, df)

        for fra in fras:
            endDate = fra._maturityDate
            df = liborCurve.df(endDate)
            print(endDate, df)

        for swap in swaps:
            endDate = swap._maturityDate
            df = liborCurve.df(endDate)
            print(endDate, df)

        swap.printFixedLegPV(spotDate)
        swap.printFloatLegPV(spotDate)

#        print(liborCurve)

###############################################################################


def test_derivativePricingExample():

    valuationDate = FinDate(10, 11, 2011)

    dccType = FinDayCountTypes.ACT_360
    depos = []

    # We do the O/N rate which settles on trade date
    spotDays = 0
    settlementDate = valuationDate.addWorkDays(spotDays)

    depositRate = 0.001410
    depo = FinLiborDeposit(settlementDate, "ON", depositRate, dccType)
    depos.append(depo)

    spotDays = 1
    settlementDate = valuationDate.addWorkDays(spotDays)

    depositRate = 0.001410
    depo = FinLiborDeposit(settlementDate, "TN", depositRate, dccType)
    depos.append(depo)

    spotDays = 2
    settlementDate = valuationDate.addWorkDays(spotDays)

    depositRate = 0.001910
    depo = FinLiborDeposit(settlementDate, "1W", depositRate, dccType)
    depos.append(depo)

    depositRate = 0.002090
    depo = FinLiborDeposit(settlementDate, "2W", depositRate, dccType)
    depos.append(depo)

    depositRate = 0.002490
    depo = FinLiborDeposit(settlementDate, "1M", depositRate, dccType)
    depos.append(depo)

    depositRate = 0.003450
    depo = FinLiborDeposit(settlementDate, "2M", depositRate, dccType)
    depos.append(depo)

    depositRate = 0.004570
    depo = FinLiborDeposit(settlementDate, "3M", depositRate, dccType)
    depos.append(depo)

    depositRate = 0.005230
    depo = FinLiborDeposit(settlementDate, "4M", depositRate, dccType)
    depos.append(depo)

    depositRate = 0.005860
    depo = FinLiborDeposit(settlementDate, "5M", depositRate, dccType)
    depos.append(depo)

    depositRate = 0.006540
    depo = FinLiborDeposit(settlementDate, "6M", depositRate, dccType)
    depos.append(depo)

    depositRate = 0.007080
    depo = FinLiborDeposit(settlementDate, "7M", depositRate, dccType)
    depos.append(depo)

    depositRate = 0.007540
    depo = FinLiborDeposit(settlementDate, "8M", depositRate, dccType)
    depos.append(depo)

    depositRate = 0.008080
    depo = FinLiborDeposit(settlementDate, "9M", depositRate, dccType)
    depos.append(depo)

    depositRate = 0.008570
    depo = FinLiborDeposit(settlementDate, "10M", depositRate, dccType)
    depos.append(depo)

    depositRate = 0.009130
    depo = FinLiborDeposit(settlementDate, "11M", depositRate, dccType)
    depos.append(depo)

    fras = []

    swaps = []
    dayCountType = FinDayCountTypes.THIRTY_E_360_ISDA
#    dayCountType = FinDayCountTypes.ACT_360
    freqType = FinFrequencyTypes.SEMI_ANNUAL

    swapRate = 0.0058
    swap = FinLiborSwap(settlementDate, "1Y", swapRate, freqType, dayCountType)
    swaps.append(swap)

    swapRate = 0.0060
    swap = FinLiborSwap(settlementDate, "2Y", swapRate, freqType, dayCountType)
    swaps.append(swap)

    swapRate = 0.0072
    swap = FinLiborSwap(settlementDate, "3Y", swapRate, freqType, dayCountType)
    swaps.append(swap)

    swapRate = 0.0096
    swap = FinLiborSwap(settlementDate, "4Y", swapRate, freqType, dayCountType)
    swaps.append(swap)

    swapRate = 0.0124
    swap = FinLiborSwap(settlementDate, "5Y", swapRate, freqType, dayCountType)
    swaps.append(swap)

    swapRate = 0.0173
    swap = FinLiborSwap(settlementDate, "7Y", swapRate, freqType, dayCountType)
    swaps.append(swap)

    swapRate = 0.0219
    swap = FinLiborSwap(settlementDate, "10Y", swapRate, freqType, dayCountType)
    swaps.append(swap)

    swapRate = 0.0283
    swap = FinLiborSwap(settlementDate, "30Y", swapRate, freqType, dayCountType)
    swaps.append(swap)

    numRepeats = 10
    start = time.time()

    for _ in range(0, numRepeats):
        _ = FinLiborCurve(valuationDate, depos, fras, swaps,
                                   FinInterpTypes.FLAT_FORWARDS)

    end = time.time()
    elapsed1 = end - start

    start = time.time()

    for _ in range(0, numRepeats):
        _ = FinLiborCurve(valuationDate, depos, fras, swaps,
                                   FinInterpTypes.LINEAR_SWAP_RATES)

    end = time.time()
    elapsed2 = end - start

    testCases.header("METHOD", "TIME")
    testCases.print("NON-LINEAR SOLVER BOOTSTRAP", elapsed1/numRepeats)
    testCases.print("LINEAR SWAP BOOTSTRAP", elapsed2/numRepeats)

###############################################################################


def test_bloombergPricingExample():

    ''' This is an example of a replication of a BBG example from
    https://github.com/vilen22/curve-building/blob/master/Bloomberg%20Curve%20Building%20Replication.xlsx

    '''
    valuationDate = FinDate(6, 6, 2018)

    # We do the O/N rate which settles on trade date
    spotDays = 0
    settlementDate = valuationDate.addWorkDays(spotDays)
    depoDCCType = FinDayCountTypes.ACT_360
    depos = []
    depositRate = 0.0231381
    maturityDate = settlementDate.addMonths(3)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate,
                           depoDCCType)
    depos.append(depo)

    futs = []
    fut = FinLiborFuture(valuationDate, 1); futs.append(fut)
    fut = FinLiborFuture(valuationDate, 2); futs.append(fut)
    fut = FinLiborFuture(valuationDate, 3); futs.append(fut)
    fut = FinLiborFuture(valuationDate, 4); futs.append(fut)
    fut = FinLiborFuture(valuationDate, 5); futs.append(fut)
    fut = FinLiborFuture(valuationDate, 6); futs.append(fut)

    fras = [None]*6
    fras[0] = futs[0].toFRA(97.6675, -0.00005)
    fras[1] = futs[1].toFRA(97.5200, -0.00060)
    fras[2] = futs[2].toFRA(97.3550, -0.00146)
    fras[3] = futs[3].toFRA(97.2450, -0.00263)
    fras[4] = futs[4].toFRA(97.1450, -0.00411)
    fras[5] = futs[5].toFRA(97.0750, -0.00589)

    accrual = FinDayCountTypes.THIRTY_E_360
    freq = FinFrequencyTypes.SEMI_ANNUAL

    spotDays = 2
    settlementDate = valuationDate.addWorkDays(spotDays)
    notional = ONE_MILLION

    swaps = []
    swap = FinLiborSwap(settlementDate, "2Y", (2.77417+2.77844)/200, freq, accrual, notional); swaps.append(swap)
    swap = FinLiborSwap(settlementDate, "3Y", (2.86098+2.86582)/200, freq, accrual); swaps.append(swap)
    swap = FinLiborSwap(settlementDate, "4Y", (2.90240+2.90620)/200, freq, accrual); swaps.append(swap)
    swap = FinLiborSwap(settlementDate, "5Y", (2.92944+2.92906)/200, freq, accrual); swaps.append(swap)
    swap = FinLiborSwap(settlementDate, "6Y", (2.94001+2.94499)/200, freq, accrual); swaps.append(swap)
    swap = FinLiborSwap(settlementDate, "7Y", (2.95352+2.95998)/200, freq, accrual); swaps.append(swap)
    swap = FinLiborSwap(settlementDate, "8Y", (2.96830+2.97400)/200, freq, accrual); swaps.append(swap)
    swap = FinLiborSwap(settlementDate, "9Y", (2.98403+2.98817)/200, freq, accrual); swaps.append(swap)
    swap = FinLiborSwap(settlementDate, "10Y", (2.99716+3.00394)/200, freq, accrual); swaps.append(swap)
    swap = FinLiborSwap(settlementDate, "11Y", (3.01344+3.01596)/200, freq, accrual); swaps.append(swap)
    swap = FinLiborSwap(settlementDate, "12Y", (3.02276+3.02684)/200, freq, accrual); swaps.append(swap)
    swap = FinLiborSwap(settlementDate, "15Y", (3.04092+3.04508)/200, freq, accrual); swaps.append(swap)
    swap = FinLiborSwap(settlementDate, "20Y", (3.04417+3.05183)/200, freq, accrual); swaps.append(swap)
    swap = FinLiborSwap(settlementDate, "25Y", (3.03219+3.03621)/200, freq, accrual); swaps.append(swap)
    swap = FinLiborSwap(settlementDate, "30Y", (3.01030+3.01370)/200, freq, accrual); swaps.append(swap)
    swap = FinLiborSwap(settlementDate, "40Y", (2.96946+2.97354)/200, freq, accrual); swaps.append(swap)
    swap = FinLiborSwap(settlementDate, "50Y", (2.91552+2.93748)/200, freq, accrual); swaps.append(swap)

    liborCurve = FinLiborCurve(valuationDate, depos, fras, swaps)

    # The valuation of 53714.55 is very close to the spreadsheet value 53713.96
    principal = 0.0

    testCases.header("VALUATION TO TODAY DATE"," PV")
    testCases.print("VALUE:", swaps[0].value(valuationDate, liborCurve, liborCurve, None, principal))
    testCases.print("FIXED:", swaps[0].fixedLegValue(valuationDate, liborCurve, principal))
    testCases.print("FLOAT:", swaps[0].floatLegValue(valuationDate, liborCurve, liborCurve, None, principal))

    testCases.header("VALUATION TO SWAP SETTLEMENT DATE"," PV")
    testCases.print("VALUE:", swaps[0].value(settlementDate, liborCurve, liborCurve, None, principal))
    testCases.print("FIXED:", swaps[0].fixedLegValue(settlementDate, liborCurve, principal))
    testCases.print("FLOAT:", swaps[0].floatLegValue(settlementDate, liborCurve, liborCurve, None, principal))

    # swaps[0].printFixedLegPV()
    # swaps[0].printFloatLegPV()

###############################################################################


test_bloombergPricingExample()
test_derivativePricingExample()
test_FinLiborDepositsOnly()
test_FinLiborFRAsOnly()
test_FinLiborDepositsFRAsSwaps()
test_FinLiborDepositsFuturesSwaps()

testCases.compareTestCases()
