# -*- coding: utf-8 -*-

import numpy as np
from financepy.finutils.FinDate import FinDate
from financepy.models.FinHullWhiteRateModel import FinHullWhiteRateModel
from financepy.models.FinHullWhiteRateModel import printTree
from financepy.market.curves.FinDiscountCurve import FinDiscountCurve
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
    vAnal = model.optionOnZeroCouponBond_Anal(startDate,
                                              expiryDate,
                                              maturityDate,
                                              strike, face,
                                              curve)

    # Test convergence
    numStepsList = range(50, 2000, 25)
    analVector = []
    treeVector = []

    for numTimeSteps in numStepsList:
        start = time.time()
        model.buildTree(startDate, expiryDate, numTimeSteps, curve)
        vTree = model.optionOnZeroCouponBond_Tree(startDate, expiryDate,
                                                  maturityDate, strike, face)
        end = time.time()
        period = end-start
        treeVector.append(vTree[1])
        analVector.append(vAnal[1])
        print(numTimeSteps, vTree, vAnal, period)

    plt.plot(numStepsList, treeVector)
    plt.plot(numStepsList, analVector)

###############################################################################


test_HullWhiteExampleOne()
test_HullWhiteExampleTwo()
