###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
import matplotlib.pyplot as plt
import numpy as np

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinCalendar import FinCalendarTypes
from financepy.products.libor.FinLiborCurve import FinLiborCurve
from financepy.products.libor.FinLiborFRA import FinLiborFRA
from financepy.products.libor.FinLiborDeposit import FinLiborDeposit
from financepy.products.libor.FinLiborSwap import FinLiborSwap
from financepy.finutils.FinCalendar import FinBusDayAdjustTypes

sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)


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

    liborCurve = FinLiborCurve("USD_LIBOR",
                               settlementDate,
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

    liborCurve = FinLiborCurve("USD_LIBOR",
                               settlementDate,
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
    fixedDCCType = FinDayCountTypes.ACT_365_ISDA
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

    liborCurve = FinLiborCurve("USD_LIBOR",
                               valuationDate,
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

    # times = np.linspace(0, 10, 100)
    # fwds = liborCurve.fwd(times)
    # plt.plot(times, fwds)

    if 1 == 0:
        for depo in depos:
            v = depo.value(valuationDate, liborCurve)
            df = liborCurve.df(depo._maturityDate)
            print("Depo Value:", v, df)

        for fra in fras:
            v = fra.value(valuationDate, liborCurve)
            df = liborCurve.df(fra._maturityDate)
            print("FRA Value:", v, df)

        for swap in swaps:
            v = swap.value(valuationDate, liborCurve, liborCurve)
            df = liborCurve.df(swap._maturityDate)
            print(swap._maturityDate, "Swap Value:", v, df)

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
    fixedDCCType = FinDayCountTypes.THIRTY_360
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

    liborCurve = FinLiborCurve("USD_LIBOR", spotDate, depos, fras, swaps)

    times = np.linspace(0.0, 2.0, 25)
    dates = spotDate.addYears(times)
    zeroRates = liborCurve.zeroRate(dates)
    fwdRates = liborCurve.fwd(dates)

    if 1 == 0:
        plt.figure(figsize=(8, 6))
        plt.plot(times, zeroRates*100, label="zero rates")
        plt.plot(times, fwdRates*100, label="fwd rates")
        plt.xlabel("Times")
        plt.ylabel("CC forward rates")
        plt.legend()

        print("==============================================================")
        for fra in fras:
            fra.print()
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

        swap.printFixedLeg(spotDate)
        swap.printFloatLeg(spotDate)

        liborCurve.print()

###############################################################################


test_FinLiborDepositsOnly()
test_FinLiborFRAsOnly() 
test_FinLiborDepositsFRAsSwaps()
test_FinLiborDepositsFuturesSwaps()

testCases.compareTestCases()
