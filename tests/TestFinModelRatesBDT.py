##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("..")

from financepy.finutils.FinDate import FinDate
from financepy.market.curves.FinDiscountCurve import FinDiscountCurve
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.products.bonds.FinBond import FinBond
from financepy.products.rates.FinIborSwaption import FinIborSwaption
from financepy.products.rates.FinIborSwaption import FinSwapTypes
from financepy.models.FinModelBlack import FinModelBlack
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinGlobalVariables import gDaysInYear
from financepy.market.curves.FinDiscountCurveZeros import FinDiscountCurveZeros
from financepy.models.FinModelRatesBDT import FinModelRatesBDT
from financepy.finutils.FinHelperFunctions import printTree
from financepy.finutils.FinGlobalTypes import FinExerciseTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

PLOT_GRAPHS = False

###############################################################################

def testBlackModelCheck():

    # Checking Andersen paper using Black's model
    # Used to check swaption price below - we have Ts = 1 and Te = 4
    # Expect a price around 122 cents which is what I find.

    valuationDate = FinDate(1, 1, 2020)
    liborCurve = FinDiscountCurveFlat(valuationDate, 0.06,
                                      FinFrequencyTypes.SEMI_ANNUAL)

    settlementDate = FinDate(1, 1, 2020)
    exerciseDate = FinDate(1, 1, 2021)
    maturityDate = FinDate(1, 1, 2024)

    fixedCoupon = 0.06
    fixedFrequencyType = FinFrequencyTypes.SEMI_ANNUAL
    fixedDayCountType = FinDayCountTypes.THIRTY_E_360_ISDA
    notional = 100.0

    # Pricing a PAY
    swaptionType = FinSwapTypes.PAY
    swaption = FinIborSwaption(settlementDate,
                                exerciseDate,
                                maturityDate,
                                swaptionType,
                                fixedCoupon,
                                fixedFrequencyType,
                                fixedDayCountType,
                                notional)

    model = FinModelBlack(0.20)
    v = swaption.value(valuationDate, liborCurve, model)
    testCases.header("LABEL", "VALUE")
    testCases.print("BLACK'S MODEL PRICE:", v*100)

###############################################################################


def test_BDTExampleOne():
    # HULL BOOK NOTES
    # http://www-2.rotman.utoronto.ca/~hull/technicalnotes/TechnicalNote23.pdf

    valuationDate = FinDate(1, 1, 2020)
    years = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    zeroDates = valuationDate.addYears(years)
    zeroRates = [0.00, 0.10, 0.11, 0.12, 0.125, 0.13]

    testCases.header("DATES")
    testCases.print(zeroDates)

    testCases.header("RATES")
    testCases.print(zeroRates)

    curve = FinDiscountCurveZeros(valuationDate,
                                  zeroDates,
                                  zeroRates,
                                  FinFrequencyTypes.ANNUAL)

    yieldVol = 0.16

    numTimeSteps = 5
    tmat = years[-1]
    dfs = curve.df(zeroDates)

    testCases.print("DFS")
    testCases.print(dfs)

    years = np.array(years)
    dfs = np.array(dfs)

    model = FinModelRatesBDT(yieldVol, numTimeSteps)
    model.buildTree(tmat, years, dfs)

###############################################################################


def test_BDTExampleTwo():
    # Valuation of a European option on a coupon bearing bond
    # This follows example in Fig 28.11 of John Hull's book (6th Edition)
    # but does not have the exact same dt so there are some differences

    testCases.banner("===================== FIG 28.11 HULL BOOK =============")

    settlementDate = FinDate(1, 12, 2019)
    issueDate = FinDate(1, 12, 2015)
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

    strikePrice = 105.0
    face = 100.0

    tmat = (maturityDate - settlementDate) / gDaysInYear
    texp = (expiryDate - settlementDate) / gDaysInYear
    times = np.linspace(0, tmat, 11)
    dates = settlementDate.addYears(times)
    dfs = np.exp(-0.05*times)

    testCases.header("LABEL", "VALUES")
    testCases.print("TIMES:", times)

    curve = FinDiscountCurve(settlementDate, dates, dfs)

    price = bond.cleanPriceFromDiscountCurve(settlementDate, curve)
    testCases.print("Fixed Income Price:", price)

    sigma = 0.20

    # Test convergence
    numStepsList = [5] #[100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    exerciseType = FinExerciseTypes.AMERICAN

    testCases.header("Values")
    treeVector = []
    for numTimeSteps in numStepsList:
        model = FinModelRatesBDT(sigma, numTimeSteps)
        model.buildTree(tmat, times, dfs)
        v = model.bondOption(texp, strikePrice,
                             face, couponTimes, couponFlows, exerciseType)

        testCases.print(v)
        treeVector.append(v['call'])

    if PLOT_GRAPHS:
        plt.plot(numStepsList, treeVector)

    # The value in Hull converges to 0.699 with 100 time steps while I get 0.70

    if 1 == 0:
        print("RT")
        printTree(model._rt, 5)
        print("Q")
        printTree(model._Q, 5)

###############################################################################


def test_BDTExampleThree():
    # Valuation of a swaption as in Leif Andersen's paper - see Table 1 on
    # SSRN-id155208.pdf

    testCases.banner("===================== ANDERSEN PAPER ==============")

    # This is a sanity check
    testBlackModelCheck()

    settlementDate = FinDate(1, 1, 2020)
    times = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    dates = settlementDate.addYears(times)
    rate = 0.06
    dfs = 1.0 / (1.0 + rate/2.0)**(2.0*times)
    curve = FinDiscountCurve(settlementDate, dates, dfs)

    coupon = 0.06
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    strikePrice = 100.0
    face = 100.0
    # Andersen paper
    numTimeSteps = 200

    testCases.header("ExerciseType", "Sigma", "NumSteps", "Texp", "Tmat", 
                     "V_Fixed", "V_pay", "V_rec")

    for exerciseType in [FinExerciseTypes.EUROPEAN,
                         FinExerciseTypes.BERMUDAN]:

        for maturityYears in [4.0, 5.0, 10.0, 20.0]:

            maturityDate = settlementDate.addYears(maturityYears)
            issueDate = FinDate(maturityDate._d, maturityDate._m, 2000)

            if maturityYears == 4.0 or maturityYears == 5.0:
                sigma = 0.2012
            elif maturityYears == 10.0:
                sigma = 0.1522
            elif maturityYears == 20.0:
                sigma = 0.1035

            for expiryYears in range(int(maturityYears/2)-1, int(maturityYears)):

                expiryDate = settlementDate.addYears(expiryYears)

                tmat = (maturityDate - settlementDate) / gDaysInYear
                texp = (expiryDate - settlementDate) / gDaysInYear

                bond = FinBond(issueDate, maturityDate, coupon, freqType, accrualType)

                couponTimes = []
                couponFlows = []
                cpn = bond._coupon/bond._frequency
                for flowDate in bond._flowDates:
                    if flowDate > expiryDate:
                        flowTime = (flowDate - settlementDate) / gDaysInYear
                        couponTimes.append(flowTime)
                        couponFlows.append(cpn)

                couponTimes = np.array(couponTimes)
                couponFlows = np.array(couponFlows)

                price = bond.cleanPriceFromDiscountCurve(settlementDate, curve)

                model = FinModelRatesBDT(sigma, numTimeSteps)
                model.buildTree(tmat, times, dfs)

                v = model.bermudanSwaption(texp,
                                           tmat,
                                           strikePrice,
                                           face,
                                           couponTimes,
                                           couponFlows,
                                           exerciseType)

                testCases.print("%s" % exerciseType,
                                "%9.5f" % sigma,
                                "%9.5f" % numTimeSteps,
                                "%9.5f" % expiryYears,
                                "%9.5f" % maturityYears,
                                "%9.5f" % price,
                                "%9.2f" % (v['pay']*100.0),
                                "%9.2f" % (v['rec']*100.0))

###############################################################################
# This has broken and needs to be repaired!!!!


test_BDTExampleOne()
test_BDTExampleTwo()
test_BDTExampleThree()

testCases.compareTestCases()
