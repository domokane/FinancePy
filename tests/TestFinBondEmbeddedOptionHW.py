###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import matplotlib.pyplot as plt
import time

import sys
sys.path.append("..")

from financepy.utils.FinGlobalTypes import FinSwapTypes
from financepy.utils.Date import Date
from financepy.utils.Frequency import FinFrequencyTypes
from financepy.utils.DayCount import FinDayCountTypes
from financepy.models.FinModelRatesHW import FinModelRatesHW

from financepy.products.rates.IborSwap import FinIborSwap
from financepy.products.rates.FinIborSingleCurve import FinIborSingleCurve
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.products.bonds.FinBond import FinBond
from financepy.products.bonds.FinBondEmbeddedOption import FinBondEmbeddedOption

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

plotGraphs = False

###############################################################################


def test_FinBondEmbeddedOptionMATLAB():

    # https://fr.mathworks.com/help/fininst/optembndbyhw.html
    # I FIND THAT THE PRICE CONVERGES TO 102.88 WHICH IS CLOSE TO 102.9127
    # FOUND BY MATLAB ALTHOUGH THEY DO NOT EXAMINE THE ASYMPTOTIC PRICE
    # WHICH MIGHT BE A BETTER MATCH

    settlement_date = Date(1, 1, 2007)
    valuation_date = settlement_date
    
    ###########################################################################

    dcType = FinDayCountTypes.THIRTY_E_360
    fixedFreq = FinFrequencyTypes.ANNUAL
    fixed_legType = FinSwapTypes.PAY
    swap1 = FinIborSwap(settlement_date, "1Y", fixed_legType, 0.0350, fixedFreq, dcType)
    swap2 = FinIborSwap(settlement_date, "2Y", fixed_legType, 0.0400, fixedFreq, dcType)
    swap3 = FinIborSwap(settlement_date, "3Y", fixed_legType, 0.0450, fixedFreq, dcType)
    swaps = [swap1, swap2, swap3]
    discount_curve = FinIborSingleCurve(valuation_date, [], [], swaps)

    ###########################################################################

    issue_date = Date(1, 1, 2004)
    maturity_date = Date(1, 1, 2010)
    
    coupon = 0.0525
    freq_type = FinFrequencyTypes.ANNUAL
    accrual_type = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(issue_date, maturity_date, coupon, freq_type, accrual_type)

    callDates = []
    callPrices = []
    putDates = []
    putPrices = []

    putDate = Date(1, 1, 2008)
    for _ in range(0, 24):
        putDates.append(putDate)
        putPrices.append(100)
        putDate = putDate.addMonths(1)

    testCases.header("BOND PRICE", "PRICE")
    v = bond.cleanPriceFromDiscountCurve(settlement_date, discount_curve)
    testCases.print("Bond Pure Price:", v)

    sigma = 0.01  # basis point volatility
    a = 0.1

    puttableBond = FinBondEmbeddedOption(issue_date,
                                         maturity_date, coupon,
                                         freq_type, accrual_type,
                                         callDates, callPrices,
                                         putDates, putPrices)

    testCases.header("TIME", "NumTimeSteps", "BondWithOption", "BondPure")

    timeSteps = range(50, 1000, 10)
    values = []
    for numTimeSteps in timeSteps:
        model = FinModelRatesHW(sigma, a, numTimeSteps)
        start = time.time()
        v = puttableBond.value(settlement_date, discount_curve, model)
        end = time.time()
        period = end - start
        testCases.print(period, numTimeSteps, v['bondwithoption'],
                        v['bondpure'])
        values.append(v['bondwithoption'])

    if plotGraphs:
        plt.figure()
        plt.plot(timeSteps, values)

###############################################################################


def test_FinBondEmbeddedOptionQUANTLIB():

    # Based on example at the nice blog on Quantlib at
    # http://gouthamanbalaraman.com/blog/callable-bond-quantlib-python.html
    # I get a price of 68.97 for 1000 time steps which is higher than the
    # 68.38 found in blog article. But this is for 40 grid points.
    # Note also that a basis point vol of 0.120 is 12% which is VERY HIGH!

    valuation_date = Date(16, 8, 2016)
    settlement_date = valuation_date.addWeekDays(3)

    ###########################################################################

    discount_curve = FinDiscountCurveFlat(valuation_date, 0.035,
                                 FinFrequencyTypes.SEMI_ANNUAL)

    ###########################################################################

    issue_date = Date(15, 9, 2010)
    maturity_date = Date(15, 9, 2022)
    coupon = 0.025
    freq_type = FinFrequencyTypes.QUARTERLY
    accrual_type = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(issue_date, maturity_date, coupon, freq_type, accrual_type)

    ###########################################################################
    # Set up the call and put times and prices
    ###########################################################################

    nextCallDate = Date(15, 9, 2016)
    callDates = [nextCallDate]
    callPrices = [100.0]

    for _ in range(1, 24):
        nextCallDate = nextCallDate.addMonths(3)
        callDates.append(nextCallDate)
        callPrices.append(100.0)

    putDates = []
    putPrices = []

    # the value used in blog of 12% bp vol is unrealistic
    sigma = 0.12  # basis point volatility
    a = 0.03

    puttableBond = FinBondEmbeddedOption(issue_date,
                                         maturity_date, coupon,
                                         freq_type, accrual_type,
                                         callDates, callPrices,
                                         putDates, putPrices)

    testCases.header("BOND PRICE", "PRICE")
    v = bond.cleanPriceFromDiscountCurve(settlement_date, discount_curve)
    testCases.print("Bond Pure Price:", v)

    testCases.header("TIME", "NumTimeSteps", "BondWithOption", "BondPure")
    timeSteps = range(100, 1000, 100)
    values = []
    for numTimeSteps in timeSteps:
        model = FinModelRatesHW(sigma, a, numTimeSteps)
        start = time.time()
        v = puttableBond.value(settlement_date, discount_curve, model)
        end = time.time()
        period = end - start
        testCases.print(period, numTimeSteps, v['bondwithoption'], v['bondpure'])
        values.append(v['bondwithoption'])

    if plotGraphs:
        plt.figure()
        plt.title("Puttable Bond Price Convergence")
        plt.plot(timeSteps, values)

###############################################################################


test_FinBondEmbeddedOptionMATLAB()
test_FinBondEmbeddedOptionQUANTLIB()
testCases.compareTestCases()
