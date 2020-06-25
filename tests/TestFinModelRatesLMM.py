# -*- coding: utf-8 -*-

import numpy as np
import time

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.finutils.FinDate import FinDate
from financepy.models.FinModelRatesLMM import FinModelRatesLMM
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.market.curves.FinZeroCurve import FinZeroCurve
from financepy.market.curves.FinInterpolate import FinInterpMethods
from financepy.finutils.FinCalendar import FinCalendarTypes
from financepy.finutils.FinCalendar import FinBusDayAdjustTypes
from financepy.finutils.FinCalendar import FinDateGenRuleTypes
from financepy.products.libor.FinLiborCapFloor import FinLiborCapFloorTypes
from financepy.products.libor.FinLiborCapFloor import FinLiborCapFloor
from financepy.models.FinModelBlack import FinModelBlack
from financepy.models.FinModelRatesLMM import FinModelRatesLMM
from financepy.finutils.FinSchedule import FinSchedule
from financepy.market.volatility.FinCapletVolCurve import FinCapletVolCurve

testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_LMMLiborsOne():

    valuationDate = FinDate(20, 6, 2016)

    numMaturities = 40
    times = np.linspace(0.0, 10.0, numMaturities+1)
    zeros = np.linspace(0.05, 0.05, numMaturities+1)
    liborCurve = FinZeroCurve(valuationDate, times, zeros)

    numMaturities = len(times)

    start = time.time()

    startDate = valuationDate
    endDate = startDate.addTenor("2Y")
    freqType = FinFrequencyTypes.ANNUAL
    dayCountType = FinDayCountTypes.ACT_ACT_ISDA
    schedule = FinSchedule(valuationDate, endDate, freqType)
    paymentDates = schedule._adjustedDates

    lmm = FinModelRatesLMM(valuationDate,
                           paymentDates,
                           dayCountType)

    capletDates = paymentDates[1:]

    vol = 0.2
    capletVols = np.linspace(vol, vol, len(capletDates))

    print(capletDates)
    print(capletVols)

    capletVolCurve = FinCapletVolCurve(valuationDate,
                                       capletDates,
                                       capletVols,
                                       dayCountType)

    numPaths = 10000
    seed = 4278

    libors = lmm.simulateLibors(liborCurve, capletVolCurve,
                                numPaths, seed)
    end = time.time()
    period = end - start
    print("Time taken:", period)

    for t in lmm._fwdEndTimes:
        df1 = lmm.bondPrice(t)
        df2 = liborCurve.df(t)
        df3 = np.exp(-t*0.05)
        print("t:", t, "df LMM:", df1, "df CURVE:", df2, "df FLAT:", df3)

#    for iPath in range(0, numPaths):
#        for iTime in range(0, len(sigmas)):
#            print("Path:",iPath,"Time:",iTime,"LIBORS:",libors[:,iTime,iPath])

###############################################################################


def test_LiborMarketModelExampleTwo():

    valuationDate = FinDate(14, 6, 2016)

    dates = [FinDate(14, 6, 2016), FinDate(14, 9, 2016),
             FinDate(14, 12, 2016), FinDate(14, 6, 2017),
             FinDate(14, 6, 2019), FinDate(14, 6, 2021),
             FinDate(15, 6, 2026), FinDate(16, 6, 2031),
             FinDate(16, 6, 2036), FinDate(14, 6, 2046)]

    rates = [0.000000, 0.006616, 0.007049, 0.007795,
             0.009599, 0.011203, 0.015068, 0.017583,
             0.018998, 0.020080]

    frequencyType = FinFrequencyTypes.ANNUAL
    dayCountType = FinDayCountTypes.ACT_ACT_ISDA

    discountCurve = FinZeroCurve(valuationDate, dates, rates, frequencyType,
                                 dayCountType,
                                 FinInterpMethods.LINEAR_ZERO_RATES)

    startDate = FinDate(14, 6, 2016)
    endDate = FinDate(14, 6, 2026)
    calendarType = FinCalendarTypes.US
    busDayAdjustType = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    dateGenRuleType = FinDateGenRuleTypes.FORWARD
    lastFixing = 0.0065560
    notional = 1000000
    dayCountType = FinDayCountTypes.ACT_360
    optionType = FinLiborCapFloorTypes.CAP
    strikeRate = 0.02

    cap = FinLiborCapFloor(startDate, endDate, optionType, strikeRate,
                           lastFixing, frequencyType,  dayCountType, notional,
                           calendarType, busDayAdjustType, dateGenRuleType)

    blackVol = 0.547295
    modelBlack = FinModelBlack(blackVol)

    start = time.time()
    v1 = cap.value(valuationDate, discountCurve, modelBlack)
    end = time.time()
    period = end - start
    print("V1:", v1, period)

    sigma = blackVol
    paymentDates = cap._capFloorDates
    numPaths = 10
    seed = 42

    start = time.time()

    modelLMM = FinModelRatesLMM(valuationDate)
    modelLMM.simulateLibors(paymentDates,
                            discountCurve,
                            dayCountType,
                            sigma,
                            numPaths,
                            seed)

    v2 = modelLMM.valueCapFloor(strikeRate, optionType, notional)
    end = time.time()
    period = end - start

    print("V2:", v2, period)

###############################################################################

def test_LMMLiborsHullExample():

    valuationDate = FinDate(1, 1, 2020)
    numMaturities = 10
    times = np.linspace(0.0, 10.0, numMaturities+1)
    zeros = np.linspace(0.05, 0.05, numMaturities+1)

    numMaturities = len(times)
    liborCurve = FinZeroCurve(valuationDate, times, zeros)

    numPaths = 10000
    seed = 42

    start = time.time()

    startDate = valuationDate
    endDate = startDate.addTenor("10Y")
    freqType = FinFrequencyTypes.ANNUAL
    dayCountType = FinDayCountTypes.ACT_ACT_ISDA
    schedule = FinSchedule(valuationDate, endDate, freqType)
    paymentDates = schedule._adjustedDates

    capletDates = paymentDates[1:]
#    capletVols = np.array([15.50, 18.25, 17.91, 17.74, 17.27,
#                  16.79, 16.30, 16.01, 15.76, 15.54]) / 100.0

    vol = 0.40
    capletVols = np.linspace(vol, vol, 10)

    capletVolCurve = FinCapletVolCurve(valuationDate,
                                       capletDates,
                                       capletVols,
                                       dayCountType)

    capletVolCurve.dump()

    lmm = FinModelRatesLMM(valuationDate,
                           paymentDates,
                           dayCountType)

    libors = lmm.simulateLibors(liborCurve, capletVolCurve, numPaths, seed)

    end = time.time()
    period = end - start
    print("Time taken:", period)

#    lmm.dump()

    times = lmm._fwdEndTimes

    for t in times:
        df1 = lmm.bondPrice(t)
        df2 = liborCurve.df(t)
        df3 = np.exp(-0.05*t)
        print("t:", t, "df LMM:", df1, "df CURVE:", df2, "FLAT DF:", df3)

#    for iPath in range(0, numPaths):
#        for iTime in range(0, len(sigmas)):
#            print("Path:",iPath,"Time:",iTime,"LIBORS:",libors[:

###############################################################################

test_LMMLiborsOne()
# test_LiborMarketModelExampleTwo()
# test_LMMLiborsHullExample()
