# -*- coding: utf-8 -*-

import numpy as np
from financepy.finutils.FinDate import FinDate
from financepy.models.FinHullWhiteRateModel import FinHullWhiteRateModel
from financepy.market.curves.FinDiscountCurve import FinDiscountCurve
from financepy.products.bonds.FinBond import FinBond
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinGlobalVariables import gDaysInYear
from financepy.finutils.FinHelperFunctions import printTree

import matplotlib.pyplot as plt
import time

###############################################################################


def test_HullWhiteExampleOne():
    # HULL BOOK INITIAL EXAMPLE SECTION 28.7 HW EDITION 6

    times = [0.0, 0.5000, 1.00000, 1.50000, 2.00000, 2.500000, 3.00000]
    zeros = [0.03, 0.0343, 0.03824, 0.04183, 0.04512, 0.048512, 0.05086]
    times = np.array(times)
    zeros = np.array(zeros)
    dfs = np.exp(-zeros*times)

    startDate = FinDate(1, 12, 2019)
    curve = FinDiscountCurve(startDate, times, dfs)
    endDate = FinDate(1, 12, 2022)
    sigma = 0.01
    a = 0.1
    numTimeSteps = 3
    model = FinHullWhiteRateModel(a, sigma)
    model.buildTree(startDate, endDate, numTimeSteps, curve)
    printTree(model._Q)
    print("")
    printTree(model._rt)
    print("")

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
    curve = FinDiscountCurve(startDate, times, dfs)
    sigma = 0.01
    a = 0.1
    strike = 63.0
    face = 100.0

    expiryDate = startDate.addTenor("3Y")
    maturityDate = startDate.addTenor("9Y")

    model = FinHullWhiteRateModel(a, sigma)
    vAnal = model.europeanOptionOnZeroCouponBond_Anal(startDate,
                                                      expiryDate,
                                                      maturityDate,
                                                      strike, face,
                                                      curve)

    # Test convergence
    numStepsList = range(100, 500, 10)
    analVector = []
    treeVector = []

    for numTimeSteps in numStepsList:
        start = time.time()
        model.buildTree(startDate, expiryDate, numTimeSteps, curve)
        vTree = model.europeanOptionOnZeroCouponBond_Tree(startDate,
                                                          expiryDate,
                                                          maturityDate,
                                                          strike, face)
        end = time.time()
        period = end-start
        treeVector.append(vTree[1])
        analVector.append(vAnal[1])
        print(numTimeSteps, vTree, vAnal, period)

    plt.plot(numStepsList, treeVector)
    plt.plot(numStepsList, analVector)

###############################################################################


def test_HullWhiteExampleThree():
    # Valuation of a European option on a coupon bearing bond

    settlementDate = FinDate(1, 12, 2019)
    expiryDate = settlementDate.addTenor("18m")
    maturityDate = settlementDate.addTenor("10Y")
    coupon = 0.05
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(maturityDate, coupon, frequencyType, accrualType)

    strikePrice = 105.0
    face = 100.0

    tmat = (maturityDate - settlementDate) / gDaysInYear
    times = np.linspace(0, tmat, 20)
    dfs = np.exp(-0.05*times)
    curve = FinDiscountCurve(settlementDate, times, dfs)

    price = bond.fullPriceFromDiscountCurve(settlementDate, curve)
    print("Fixed Income Price:", price)

    sigma = 0.01
    a = 0.1

    # Test convergence
    numStepsList = [100,101,200,300,400,500,600,700,800,900,1000]
    isAmerican = True

    treeVector = []
    for numTimeSteps in numStepsList:
        start = time.time()
        model = FinHullWhiteRateModel(a, sigma)
        model.buildTree(settlementDate, expiryDate, numTimeSteps, curve)
        v = model.bondOption(settlementDate, expiryDate, strikePrice,
                             face, bond, isAmerican)
        end = time.time()
        period = end-start
        treeVector.append(v[0])
        print(numTimeSteps, v, period)

#        print("BOND")
#        printTree(model._bondValues)
#        print("OPTION")
#        printTree(model._optionValues)

    plt.plot(numStepsList, treeVector)

    if 1==0:
        print("RT")
        printTree(model._rt, 5)
        print("BOND")
        printTree(model._bondValues, 5)
        print("OPTION")
        printTree(model._optionValues, 5)

###############################################################################


#test_HullWhiteExampleOne()
#test_HullWhiteExampleTwo()
test_HullWhiteExampleThree()
