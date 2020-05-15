# -*- coding: utf-8 -*-

import numpy as np
import time
import matplotlib.pyplot as plt

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.models.FinModelRatesHW import FinModelRatesHW
from financepy.models.FinModelRatesBK import FinModelRatesBK

from financepy.products.libor.FinLiborSwap import FinLiborSwap
from financepy.market.curves.FinLiborCurve import FinLiborCurve
from financepy.market.curves.FinFlatCurve import FinFlatCurve
from financepy.products.bonds.FinBond import FinBond
from financepy.products.bonds.FinBondEmbeddedOption import FinBondEmbeddedOption

###############################################################################


def test_FinBondEmbeddedOptionMATLAB():

    # Based on example
    # https://fr.mathworks.com/help/fininst/optembndbyhw.html#bviuizn-1_sep_optembndbyhw_example1
    settlementDate = FinDate(1, 1, 2007)

    ###########################################################################

    dcType = FinDayCountTypes.ACT_360
    fixedFreq = FinFrequencyTypes.SEMI_ANNUAL
    swap1 = FinLiborSwap(settlementDate, FinDate(1,1,2008), 0.05, fixedFreq, dcType)
    swap2 = FinLiborSwap(settlementDate, FinDate(1,1,2009), 0.05, fixedFreq, dcType)
    swap3 = FinLiborSwap(settlementDate, FinDate(1,1,2010), 0.05, fixedFreq, dcType)
    swaps = [swap1, swap2, swap3]
    discountCurve = FinLiborCurve("USD_LIBOR", settlementDate, [], [], swaps)

    ###########################################################################

    maturityDate = FinDate(1, 1, 2010)
    coupon = 0.05
    frequencyType = FinFrequencyTypes.ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(maturityDate, coupon, frequencyType, accrualType)

    ###########################################################################
    # Set up the call and put times and prices
    ###########################################################################

    callDates = []
    callPrices = []
    putDates = [FinDate(1, 1, 2008), FinDate(1, 1, 2009), FinDate(1, 1, 2010)]
    putPrices = [100.0, 100.0, 100.0]

    putDates = []
    putPrices = []

    putDates = [FinDate(1, 1, 2009)]
    putPrices = [100.0]

    v = bond.valueBondUsingDiscountCurve(settlementDate, discountCurve)
    print("Pure Bond Price:", v)

    sigma = 0.01  # basis point volatility
    a = 0.1
    numTimeSteps = 5

    puttableBond = FinBondEmbeddedOption(maturityDate, coupon,
                                         frequencyType, accrualType,
                                         callDates, callPrices,
                                         putDates, putPrices)

    timeSteps = range(5, 9, 1)
    values = []
    for numTimeSteps in timeSteps:
        print("==============================================================")
        model = FinModelRatesHW(a, sigma, numTimeSteps)
        v = puttableBond.value(settlementDate, discountCurve, model)
        print("NUM_TIMES_STEPS", numTimeSteps)
        values.append(v[0])

    plt.figure()
    plt.plot(timeSteps, values)
    plt.ylim(99.5, 101)

###############################################################################


def test_FinBondEmbeddedOptionQUANTLIB():

    # Based on example
    # http://gouthamanbalaraman.com/blog/callable-bond-quantlib-python.html
    settlementDate = FinDate(16, 8, 2016)

    ###########################################################################

    discountCurve = FinFlatCurve(settlementDate, 0.035, 1)

    ###########################################################################

    maturityDate = FinDate(16, 9, 2022)
    coupon = 0.025
    frequencyType = FinFrequencyTypes.QUARTERLY
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(maturityDate, coupon, frequencyType, accrualType)
    bond.printFlows(settlementDate)

    ###########################################################################
    # Set up the call and put times and prices
    ###########################################################################

    nextCallDate = FinDate(15, 9, 2016)
    callDates = [nextCallDate]
    callPrices = [100.0]

    for i in range(1, 24):
        nextCallDate = nextCallDate.addMonths(3)
        callDates.append(nextCallDate)
        callPrices.append(100.0)

    putDates = []
    putPrices = []

    v = bond.valueBondUsingDiscountCurve(settlementDate, discountCurve)
    print("Pure Bond Price:", v)

    # the value used in blog of 12% bp vol is unrealistic
    sigma = 0.012  # basis point volatility
    a = 0.03
    numTimeSteps = 40

    model = FinModelRatesHW(a, sigma, numTimeSteps)
    puttableBond = FinBondEmbeddedOption(maturityDate, coupon,
                                         frequencyType, accrualType,
                                         callDates, callPrices,
                                         putDates, putPrices)

    v = puttableBond.value(settlementDate, discountCurve, model)

    print("BondEmbeddedOption Price:", v)

    timeSteps = range(20, 300, 1)
    values = []
    for numTimeSteps in timeSteps:
        model = FinModelRatesHW(a, sigma, numTimeSteps)
        v = puttableBond.value(settlementDate, discountCurve, model)
        print(numTimeSteps, v)
        values.append(v[0])

    print(timeSteps, values)
    plt.figure()
    plt.plot(timeSteps, values)
#    plt.ylim(68, 69)

###############################################################################


test_FinBondEmbeddedOptionMATLAB()
# test_FinBondEmbeddedOptionQUANTLIB()
