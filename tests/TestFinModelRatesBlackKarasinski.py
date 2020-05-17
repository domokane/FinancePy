# -*- coding: utf-8 -*-

import numpy as np
from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinHelperFunctions import printTree
from financepy.models.FinModelRatesBK import FinModelRatesBK
from financepy.market.curves.FinDiscountCurve import FinDiscountCurve
from financepy.products.bonds.FinBond import FinBond
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinGlobalVariables import gDaysInYear
from financepy.finutils.FinHelperFunctions import printTree

import matplotlib.pyplot as plt
import time


###############################################################################

def test_BlackKarasinskiExampleOne():
    # HULL BOOK INITIAL EXAMPLE SECTION 28.7 HW EDITION 6

    times = [0.0, 0.5000, 1.00000, 1.50000, 2.00000, 2.500000, 3.00000]
    zeros = [0.03, 0.0343, 0.03824, 0.04183, 0.04512, 0.048512, 0.05086]
    times = np.array(times)
    zeros = np.array(zeros)
    dfs = np.exp(-zeros*times)

    startDate = FinDate(1, 12, 2019)
    curve = FinDiscountCurve(startDate, times, dfs)
    endDate = FinDate(1, 12, 2022)
    sigma = 0.25
    a = 0.22
    numTimeSteps = 40
    tmat = (endDate - startDate)/gDaysInYear
    model = FinModelRatesBK(a, sigma, numTimeSteps)
    model.buildTree(tmat, times, dfs)
#    printTree(model._Q)
#    print("")
#    printTree(model._rt)
#    print("")

###############################################################################


def test_BlackKarasinskiExampleTwo():
    # Valuation of a European option on a coupon bearing bond
    # This follows example in Fig 28.11 of John Hull's book but does not
    # have the exact same dt so there are some differences

    settlementDate = FinDate(1, 12, 2019)
    expiryDate = settlementDate.addTenor("18m")
    maturityDate = settlementDate.addTenor("10Y")
    coupon = 0.05
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(maturityDate, coupon, frequencyType, accrualType)

    bond.calculateFlowDates(settlementDate)
    couponTimes = []
    couponFlows = []
    cpn = bond._coupon/bond._frequency
    for flowDate in bond._flowDates:
        flowTime = (flowDate - settlementDate) / gDaysInYear
        couponTimes.append(flowTime)
        couponFlows.append(cpn)

    couponTimes = np.array(couponTimes)
    couponFlows = np.array(couponFlows)

    strikePrice = 105.0
    face = 100.0

    tmat = (maturityDate - settlementDate) / gDaysInYear
    texp = (expiryDate - settlementDate) / gDaysInYear
    times = np.linspace(0, tmat, 20)
    dfs = np.exp(-0.05*times)
    curve = FinDiscountCurve(settlementDate, times, dfs)

    price = bond.valueBondUsingDiscountCurve(settlementDate, curve)
    print("Fixed Income Price:", price)

    sigma = 0.20
    a = 0.05

    # Test convergence
    numStepsList = [100] #,101,200,300,400,500,600,700,800,900,1000]
    isAmerican = True

    treeVector = []
    for numTimeSteps in numStepsList:
        start = time.time()
        model = FinModelRatesBK(a, sigma, numTimeSteps)
        model.buildTree(tmat, times, dfs)
        v = model.bondOption(texp, strikePrice,
                             face, couponTimes, couponFlows, isAmerican)
        end = time.time()
        period = end-start
        treeVector.append(v[0])
        print(numTimeSteps, v, period)

#    plt.plot(numStepsList, treeVector)

    # The value in Hill converges to 0.699 with 100 time steps while I get 0.700

if 1==0:
    print("RT")
#    printTree(model._rt, 5)
    print("BOND")
#    printTree(model._bondValues, 5)
    print("OPTION")
#    printTree(model._optionValues, 5)

###############################################################################

# This has broken and needs to be repaired!!!!
test_BlackKarasinskiExampleOne()
# test_BlackKarasinskiExampleTwo()
