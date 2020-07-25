# -*- coding: utf-8 -*-

import numpy as np

from FinTestCases import FinTestCases, globalTestCaseMode


from financepy.finutils.FinDate import FinDate
from financepy.market.curves.FinDiscountCurve import FinDiscountCurve
from financepy.products.bonds.FinBond import FinBond
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinGlobalVariables import gDaysInYear
from financepy.market.curves.FinDiscountCurveZeros import FinDiscountCurveZeros
from financepy.models.FinModelRatesBDT import FinModelRatesBDT

import matplotlib.pyplot as plt
import time

testCases = FinTestCases(__file__, globalTestCaseMode)


###############################################################################

def test_BDTExampleOne():
    # HULL BOOK NOTES
    # http://www-2.rotman.utoronto.ca/~hull/technicalnotes/TechnicalNote23.pdf

    valuationDate = FinDate(1, 1, 2020)
    years = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    zeroDates = valuationDate.addYears(years)
    zeroRates = [0.10, 0.10, 0.11, 0.12, 0.125, 0.13]

    print(zeroDates)
    print(zeroRates)

    curve = FinDiscountCurveZeros(valuationDate,
                                  zeroDates,
                                  zeroRates,
                                  FinFrequencyTypes.ANNUAL)

    yieldVol = 0.19

    print("STARTING")
    numTimeSteps = 5
    tmat = years[-1]
    dfs = curve.df(zeroDates)

    print("DFS")
    print(dfs)

    years = np.array(years)
    dfs = np.array(dfs)

    model = FinModelRatesBDT(yieldVol, numTimeSteps)
    model.buildTree(tmat, years, dfs)

###############################################################################


def test_BDTExampleTwo():
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

    # Test convergence
    numStepsList = [100] #,101,200,300,400,500,600,700,800,900,1000]
    isAmerican = True

    treeVector = []
    for numTimeSteps in numStepsList:
        start = time.time()
        model = FinModelRatesBDT(sigma, numTimeSteps)
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
test_BDTExampleOne()
test_BDTExampleTwo()

testCases.compareTestCases()
