##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
import time

import sys
sys.path.append("..")

from financepy.utils.Date import Date
from financepy.market.curves.FinDiscountCurve import FinDiscountCurve
from financepy.products.bonds.FinBond import FinBond
from financepy.utils.Frequency import FinFrequencyTypes
from financepy.utils.DayCount import FinDayCountTypes
from financepy.utils.FinGlobalVariables import gDaysInYear
from financepy.utils.FinHelperFunctions import printTree
from financepy.models.FinModelRatesBK import FinModelRatesBK
from financepy.utils.FinGlobalTypes import FinExerciseTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

def test_BKExampleOne():

    testCases.banner("=== HULL INITIAL EXAMPLE SECTION 28.7 ED 6 PG 668 ====")

    times = [0.0, 0.5000, 1.00000, 1.50000, 2.00000, 2.500000, 3.00000]
    zeros = [0.03, 0.0343, 0.03824, 0.04183, 0.04512, 0.048512, 0.05086]
    times = np.array(times)
    zeros = np.array(zeros)
    dfs = np.exp(-zeros*times)

    start_date = Date(1, 12, 2019)
    end_date = Date(1, 6, 2021)
    sigma = 0.25
    a = 0.22
    numTimeSteps = 3
    tmat = (end_date - start_date)/gDaysInYear
    model = FinModelRatesBK(sigma, a, numTimeSteps)
    model.buildTree(tmat, times, dfs)

    # Agrees with Figure 28.10 - Not exact as we have dt not exactly 0.50
    if numTimeSteps < 5:
        testCases.header("LABEL", "VALUE")
        testCases.print("QTREE", model._Q)
        testCases.print("RTREE", model._rt)
#        printTree(model._rt)
        testCases.print("PU AT LAST TIME", model._pu)
        testCases.print("PDM AT LAST TIME", model._pm)
        testCases.print("PD AT LAST TIME", model._pd)

###############################################################################


def test_BKExampleTwo():
    # Valuation of a European option on a coupon bearing bond
    # This follows example in Fig 28.11 of John Hull's book but does not
    # have the exact same dt so there are some differences

    settlement_date = Date(1, 12, 2019)
    issue_date = Date(1, 12, 2018)
    expiry_date = settlement_date.addTenor("18m")
    maturity_date = settlement_date.addTenor("10Y")
    coupon = 0.05
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    accrual_type = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(issue_date, maturity_date, coupon, freq_type, accrual_type)

    couponTimes = []
    couponFlows = []
    cpn = bond._coupon/bond._frequency
    numFlows = len(bond._flow_dates)

    for i in range(1, numFlows):
        pcd = bond._flow_dates[i-1]
        ncd = bond._flow_dates[i]
        if pcd < settlement_date and ncd > settlement_date:
            flowTime = (pcd - settlement_date) / gDaysInYear
            couponTimes.append(flowTime)
            couponFlows.append(cpn)

    for flowDate in bond._flow_dates:
        if flowDate > settlement_date:
            flowTime = (flowDate - settlement_date) / gDaysInYear
            couponTimes.append(flowTime)
            couponFlows.append(cpn)

    couponTimes = np.array(couponTimes)
    couponFlows = np.array(couponFlows)

    strikePrice = 105.0
    face = 100.0

    tmat = (maturity_date - settlement_date) / gDaysInYear
    texp = (expiry_date - settlement_date) / gDaysInYear
    times = np.linspace(0, tmat, 11)
    dates = settlement_date.addYears(times)
    dfs = np.exp(-0.05*times)
    curve = FinDiscountCurve(settlement_date, dates, dfs)

    price = bond.cleanPriceFromDiscountCurve(settlement_date, curve)
    testCases.header("LABEL", "VALUE")
    testCases.print("Fixed Income Price:", price)

    sigma = 0.20
    a = 0.05
    numTimeSteps = 26

    model = FinModelRatesBK(sigma, a, numTimeSteps)
    model.buildTree(tmat, times, dfs)
    exerciseType = FinExerciseTypes.AMERICAN
    v = model.bondOption(texp, strikePrice, face, couponTimes,
                         couponFlows, exerciseType)

    # Test convergence
    numStepsList = [100, 200, 300, 500, 1000]
    exerciseType = FinExerciseTypes.AMERICAN

    testCases.header("TIMESTEPS", "TIME", "VALUE")
    treeVector = []
    for numTimeSteps in numStepsList:
        start = time.time()
        model = FinModelRatesBK(sigma, a, numTimeSteps)
        model.buildTree(tmat, times, dfs)
        v = model.bondOption(texp, strikePrice,
                             face, couponTimes, couponFlows, exerciseType)
        end = time.time()
        period = end-start
        treeVector.append(v)
        testCases.print(numTimeSteps, period, v)

#    plt.plot(numStepsList, treeVector)

    # Value in Hill converges to 0.699 with 100 time steps while I get 0.700

    if 1 == 0:
        print("RT")
        printTree(model._rt, 5)
        print("Q")
        printTree(model._Q, 5)

###############################################################################


test_BKExampleOne()
test_BKExampleTwo()
testCases.compareTestCases()
