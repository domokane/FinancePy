##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
import matplotlib.pyplot as plt
import time

import sys
sys.path.append("..")

from financepy.finutils.FinDate import FinDate
from financepy.models.FinModelRatesHW import FinModelRatesHW, FinHWEuropeanCalcType
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.products.bonds.FinBond import FinBond
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinGlobalVariables import gDaysInYear
from financepy.finutils.FinHelperFunctions import printTree
from financepy.finutils.FinGlobalTypes import FinExerciseTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_HullWhiteExampleOne():
    # HULL BOOK INITIAL EXAMPLE SECTION 28.7 HW EDITION 6

    times = [0.0, 0.5000, 1.00000, 1.50000, 2.00000, 2.500000, 3.00000]
    zeros = [0.03, 0.0343, 0.03824, 0.04183, 0.04512, 0.048512, 0.05086]
    times = np.array(times)
    zeros = np.array(zeros)
    dfs = np.exp(-zeros*times)

    startDate = FinDate(1, 12, 2019)
    endDate = FinDate(1, 12, 2022)
    sigma = 0.01
    a = 0.1
    numTimeSteps = 3
    model = FinModelRatesHW(sigma, a, numTimeSteps)
    treeMat = (endDate - startDate)/gDaysInYear
    model.buildTree(treeMat, times, dfs)
#   printTree(model._Q)
#   print("")
#   printTree(model._rt)
#   print("")

###############################################################################


def test_HullWhiteExampleTwo():
    # HULL BOOK ZERO COUPON BOND EXAMPLE 28.1 SEE TABLE 28.3
    # Replication may not be exact as I am using dates rather than times

    zeroDays = [0, 3, 31, 62, 94, 185, 367, 731, 1096, 1461,
                1826, 2194, 2558, 2922, 3287, 3653]

    zeroRates = [5.0, 5.01772, 4.98282, 4.97234, 4.96157, 4.99058, 5.09389,
                 5.79733, 6.30595, 6.73464, 6.94816, 7.08807, 7.27527,
                 7.30852, 7.39790, 7.49015]

    times = np.array(zeroDays) / 365.0
    zeros = np.array(zeroRates) / 100.0
    dfs = np.exp(-zeros*times)

    startDate = FinDate(1, 12, 2019)
    sigma = 0.01
    a = 0.1
    strike = 63.0
    face = 100.0

    expiryDate = startDate.addTenor("3Y")
    maturityDate = startDate.addTenor("9Y")

    texp = (expiryDate - startDate)/gDaysInYear
    tmat = (maturityDate - startDate)/gDaysInYear

    numTimeSteps = None
    model = FinModelRatesHW(sigma, a, numTimeSteps)
    vAnal = model.optionOnZCB(texp, tmat, strike, face, times, dfs)

    # Test convergence
    numStepsList = range(100, 500, 100)
    analVector = []
    treeVector = []

    testCases.banner("Comparing option on zero coupon bond analytical vs Tree")

    testCases.header("NUMTIMESTEP", "TIME", "VTREE_CALL", "VTREE_PUT", 
                     "VANAL CALL", "VANAL_PUT", "CALLDIFF", "PUTDIFF")

    for numTimeSteps in numStepsList:

        start = time.time()

        model = FinModelRatesHW(sigma, a, numTimeSteps)
        model.buildTree(texp, times, dfs)
        vTree1 = model.optionOnZeroCouponBond_Tree(texp, tmat, strike, face)

        model = FinModelRatesHW(sigma, a, numTimeSteps+1)
        model.buildTree(texp, times, dfs)
        vTree2 = model.optionOnZeroCouponBond_Tree(texp, tmat, strike, face)
 
        end = time.time()
        period = end-start
        treeVector.append(vTree1['put'])
        analVector.append(vAnal['put'])
        vTreeCall = (vTree1['call'] + vTree2['call'] ) / 2.0
        vTreePut = (vTree1['put'] + vTree2['put'] ) / 2.0
        diffC = vTreeCall - vAnal['call']
        diffP = vTreePut - vAnal['put']
        
        testCases.print(numTimeSteps, period, vTreeCall, vAnal['call'],
                        vTreePut, vAnal['put'], diffC, diffP)

 #   plt.plot(numStepsList, treeVector)
 #   plt.plot(numStepsList, analVector)

###############################################################################


def test_HullWhiteBondOption():
    # Valuation of a European option on a coupon bearing bond

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

        if ncd > settlementDate:
            
            if len(couponTimes) == 0:
                flowTime = (pcd - settlementDate) / gDaysInYear
                couponTimes.append(flowTime)
                couponFlows.append(cpn)
                
            flowTime = (ncd - settlementDate) / gDaysInYear
            couponTimes.append(flowTime)
            couponFlows.append(cpn)

    couponTimes = np.array(couponTimes)
    couponFlows = np.array(couponFlows)

    strikePrice = 100.0
    face = 100.0
    y = 0.05
    times = np.linspace(0, 10, 21)
    dfs = np.power(1+y/2, -times*2)

    sigma = 0.0000001
    a = 0.1
    model = FinModelRatesHW(sigma, a, None)

    #  Test convergence
    numStepsList = range(50, 500, 50)
    texp = (expiryDate - settlementDate)/gDaysInYear

    vJam = model.europeanBondOptionJamshidian(texp, strikePrice, face,
                                              couponTimes, couponFlows,
                                              times, dfs)
    
    testCases.banner("Pricing bond option on tree that goes to bond maturity and one using european bond option tree that goes to expiry.")

    testCases.header("NUMSTEPS", "TIME", "EXPIRY_ONLY", "EXPIRY_TREE", "JAMSHIDIAN")

    for numTimeSteps in numStepsList:

        start = time.time()
        model = FinModelRatesHW(sigma, a, numTimeSteps, FinHWEuropeanCalcType.EXPIRY_ONLY)
        model.buildTree(texp, times, dfs)

        exerciseType = FinExerciseTypes.EUROPEAN

        v1 = model.bondOption(texp, strikePrice, face,
                              couponTimes, couponFlows, exerciseType)

        model = FinModelRatesHW(sigma, a, numTimeSteps, FinHWEuropeanCalcType.EXPIRY_TREE)
        model.buildTree(texp, times, dfs)

        v2 = model.bondOption(texp, strikePrice, face,
                              couponTimes, couponFlows, exerciseType)

        end = time.time()
        period = end-start

        testCases.print(numTimeSteps, period, v1, v2, vJam)

#    plt.plot(numStepsList, treeVector)

    if 1 == 0:
        print("RT")
        printTree(model._rt, 5)
        print("BOND")
        printTree(model._bondValues, 5)
        print("OPTION")
        printTree(model._optionValues, 5)


###############################################################################


def test_HullWhiteCallableBond():
    # Valuation of a European option on a coupon bearing bond

    settlementDate = FinDate(1, 12, 2019)
    issueDate = FinDate(1, 12, 2018)
    maturityDate = settlementDate.addTenor("10Y")
    coupon = 0.05
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(issueDate, maturityDate, coupon, freqType, accrualType)

    couponTimes = []
    couponFlows = []
    cpn = bond._coupon/bond._frequency

    for flowDate in bond._flowDates[1:]:
        
        if flowDate > settlementDate:
            flowTime = (flowDate - settlementDate) / gDaysInYear
            couponTimes.append(flowTime)
            couponFlows.append(cpn)

    couponTimes = np.array(couponTimes)
    couponFlows = np.array(couponFlows)

    ###########################################################################
    # Set up the call and put times and prices
    ###########################################################################

    callDates = []
    callPrices = []
    callPx = 120.0
    callDates.append(settlementDate.addTenor("2Y")); callPrices.append(callPx)
    callDates.append(settlementDate.addTenor("3Y")); callPrices.append(callPx)
    callDates.append(settlementDate.addTenor("4Y")); callPrices.append(callPx)
    callDates.append(settlementDate.addTenor("5Y")); callPrices.append(callPx)
    callDates.append(settlementDate.addTenor("6Y")); callPrices.append(callPx)
    callDates.append(settlementDate.addTenor("7Y")); callPrices.append(callPx)
    callDates.append(settlementDate.addTenor("8Y")); callPrices.append(callPx)

    callTimes = []
    for dt in callDates:
        t = (dt - settlementDate) / gDaysInYear
        callTimes.append(t)

    putDates = []
    putPrices = []
    putPx = 98.0
    putDates.append(settlementDate.addTenor("2Y")); putPrices.append(putPx)
    putDates.append(settlementDate.addTenor("3Y")); putPrices.append(putPx)
    putDates.append(settlementDate.addTenor("4Y")); putPrices.append(putPx)
    putDates.append(settlementDate.addTenor("5Y")); putPrices.append(putPx)
    putDates.append(settlementDate.addTenor("6Y")); putPrices.append(putPx)
    putDates.append(settlementDate.addTenor("7Y")); putPrices.append(putPx)
    putDates.append(settlementDate.addTenor("8Y")); putPrices.append(putPx)

    putTimes = []
    for dt in putDates:
        t = (dt - settlementDate) / gDaysInYear
        putTimes.append(t)

    ###########################################################################

    tmat = (maturityDate - settlementDate) / gDaysInYear
    curve = FinDiscountCurveFlat(settlementDate, 0.05, FinFrequencyTypes.CONTINUOUS)

    dfs = []
    times = []

    for dt in bond._flowDates:
        if dt > settlementDate:
            t = (dt - settlementDate) / gDaysInYear
            df = curve.df(dt)
            times.append(t)
            dfs.append(df) 
                        
    dfs = np.array(dfs)
    times = np.array(times)

    ###########################################################################

    v1 = bond.cleanPriceFromDiscountCurve(settlementDate, curve)

    sigma = 0.02  # basis point volatility
    a = 0.01

    # Test convergence
    numStepsList = [100, 200, 500, 1000]
    tmat = (maturityDate - settlementDate)/gDaysInYear

    testCases.header("NUMSTEPS", "TIME", "BOND_ONLY", "CALLABLE_BOND")

    for numTimeSteps in numStepsList:

        start = time.time()
        model = FinModelRatesHW(sigma, a, numTimeSteps)
        model.buildTree(tmat, times, dfs)

        v2 = model.callablePuttableBond_Tree(couponTimes, couponFlows,
                                             callTimes, callPrices,
                                             putTimes, putPrices, 100.0)

        end = time.time()
        period = end-start
        testCases.print(numTimeSteps, period, v1, v2)

###############################################################################


test_HullWhiteExampleOne()
test_HullWhiteExampleTwo()
test_HullWhiteBondOption()
test_HullWhiteCallableBond()
testCases.compareTestCases()
