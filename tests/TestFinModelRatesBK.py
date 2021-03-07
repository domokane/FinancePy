##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
import time

import sys
sys.path.append("..")

from financepy.utils.date import Date
from financepy.market.discount.curve import DiscountCurve
from financepy.products.bonds.bond import Bond
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.global_vars import gDaysInYear
from financepy.utils.helpers import print_tree
from financepy.models.rates_bk_tree import FinModelRatesBK
from financepy.utils.global_types import FinExerciseTypes

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

    startDate = FinDate(1, 12, 2019)
    endDate = FinDate(1, 6, 2021)
    sigma = 0.25
    a = 0.22
    num_time_steps = 3
    tmat = (endDate - startDate)/gDaysInYear
    model = FinModelRatesBK(sigma, a, num_time_steps)
    model.buildTree(tmat, times, dfs)

    # Agrees with Figure 28.10 - Not exact as we have dt not exactly 0.50
    if num_time_steps < 5:
        testCases.header("LABEL", "VALUE")
        testCases.print("QTREE", model._Q)
        testCases.print("RTREE", model._rt)
#        print_tree(model._rt)
        testCases.print("PU AT LAST TIME", model._pu)
        testCases.print("PDM AT LAST TIME", model._pm)
        testCases.print("PD AT LAST TIME", model._pd)

###############################################################################


def test_BKExampleTwo():
    # Valuation of a European option on a coupon bearing bond
    # This follows example in Fig 28.11 of John Hull's book but does not
    # have the exact same dt so there are some differences

    settlementDate = FinDate(1, 12, 2019)
    issueDate = FinDate(1, 12, 2018)
    expiryDate = settlementDate.addTenor("18m")
    maturityDate = settlementDate.addTenor("10Y")
    coupon = 0.05
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(issueDate, maturityDate, coupon, freqType, accrualType)

    couponTimes = []
    couponFlows = []
    cpn = bond._coupon/bond._frequency
    numFlows = len(bond._flowDates)

    for i in range(1, numFlows):
        pcd = bond._flowDates[i-1]
        ncd = bond._flowDates[i]
        if pcd < settlementDate and ncd > settlementDate:
            flowTime = (pcd - settlementDate) / gDaysInYear
            couponTimes.append(flowTime)
            couponFlows.append(cpn)

    for flowDate in bond._flowDates:
        if flowDate > settlementDate:
            flowTime = (flowDate - settlementDate) / gDaysInYear
            couponTimes.append(flowTime)
            couponFlows.append(cpn)

    couponTimes = np.array(couponTimes)
    couponFlows = np.array(couponFlows)

    strike_price = 105.0
    face = 100.0

    tmat = (maturityDate - settlementDate) / gDaysInYear
    texp = (expiryDate - settlementDate) / gDaysInYear
    times = np.linspace(0, tmat, 11)
    dates = settlementDate.addYears(times)
    dfs = np.exp(-0.05*times)
    curve = FinDiscountCurve(settlementDate, dates, dfs)

    price = bond.cleanPriceFromDiscountCurve(settlementDate, curve)
    testCases.header("LABEL", "VALUE")
    testCases.print("Fixed Income Price:", price)

    sigma = 0.20
    a = 0.05
    num_time_steps = 26

    model = FinModelRatesBK(sigma, a, num_time_steps)
    model.buildTree(tmat, times, dfs)
    exerciseType = FinExerciseTypes.AMERICAN
    v = model.bondOption(texp, strike_price, face, couponTimes,
                         couponFlows, exerciseType)

    # Test convergence
    numStepsList = [100, 200, 300, 500, 1000]
    exerciseType = FinExerciseTypes.AMERICAN

    testCases.header("TIMESTEPS", "TIME", "VALUE")
    treeVector = []
    for num_time_steps in numStepsList:
        start = time.time()
        model = FinModelRatesBK(sigma, a, num_time_steps)
        model.buildTree(tmat, times, dfs)
        v = model.bondOption(texp, strike_price,
                             face, couponTimes, couponFlows, exerciseType)
        end = time.time()
        period = end-start
        treeVector.append(v)
        testCases.print(num_time_steps, period, v)

#    plt.plot(numStepsList, treeVector)

    # Value in Hill converges to 0.699 with 100 time steps while I get 0.700

    if 1 == 0:
        print("RT")
        print_tree(model._rt, 5)
        print("Q")
        print_tree(model._Q, 5)

###############################################################################


test_BKExampleOne()
test_BKExampleTwo()
testCases.compareTestCases()
