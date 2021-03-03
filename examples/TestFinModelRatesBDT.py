##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("..")

from financepy.utils.date import Date
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.FinDiscountCurveFlat import DiscountCurveFlat
from financepy.products.bonds.bond import Bond
from financepy.products.rates.FinIborSwaption import FinIborSwaption
from financepy.products.rates.FinIborSwaption import FinSwapTypes
from financepy.models.black import FinModelBlack
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.global_variables import gDaysInYear
from financepy.market.curves.FinDiscountCurveZeros import DiscountCurveZeros
from financepy.models.rates_bdt_tree import FinModelRatesBDT
from financepy.utils.helper_functions import printTree
from financepy.utils.FinGlobalTypes import FinExerciseTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

PLOT_GRAPHS = False

###############################################################################

def testBlackModelCheck():

    # Checking Andersen paper using Black's model
    # Used to check swaption price below - we have Ts = 1 and Te = 4
    # Expect a price around 122 cents which is what I find.

    valuation_date = Date(1, 1, 2020)
    libor_curve = DiscountCurveFlat(valuation_date, 0.06,
                                    FrequencyTypes.SEMI_ANNUAL)

    settlement_date = Date(1, 1, 2020)
    exerciseDate = Date(1, 1, 2021)
    maturity_date = Date(1, 1, 2024)

    fixedCoupon = 0.06
    fixedFrequencyType = FrequencyTypes.SEMI_ANNUAL
    fixedDayCountType = DayCountTypes.THIRTY_E_360_ISDA
    notional = 100.0

    # Pricing a PAY
    swaptionType = FinSwapTypes.PAY
    swaption = FinIborSwaption(settlement_date,
                                exerciseDate,
                                maturity_date,
                                swaptionType,
                                fixedCoupon,
                                fixedFrequencyType,
                                fixedDayCountType,
                                notional)

    model = FinModelBlack(0.20)
    v = swaption.value(valuation_date, libor_curve, model)
    testCases.header("LABEL", "VALUE")
    testCases.print("BLACK'S MODEL PRICE:", v*100)

###############################################################################


def test_BDTExampleOne():
    # HULL BOOK NOTES
    # http://www-2.rotman.utoronto.ca/~hull/technicalnotes/TechnicalNote23.pdf

    valuation_date = Date(1, 1, 2020)
    years = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    zeroDates = valuation_date.addYears(years)
    zeroRates = [0.00, 0.10, 0.11, 0.12, 0.125, 0.13]

    testCases.header("DATES")
    testCases.print(zeroDates)

    testCases.header("RATES")
    testCases.print(zeroRates)

    curve = DiscountCurveZeros(valuation_date,
                               zeroDates,
                               zeroRates,
                               FrequencyTypes.ANNUAL)

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

    settlement_date = Date(1, 12, 2019)
    issue_date = Date(1, 12, 2015)
    expiry_date = settlement_date.addTenor("18m")
    maturity_date = settlement_date.addTenor("10Y")
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_date, maturity_date, coupon, freq_type, accrual_type)

    coupon_times = []
    coupon_flows = []
    cpn = bond._coupon/bond._frequency
    num_flows = len(bond._flow_dates)

    for i in range(1, num_flows):
        pcd = bond._flow_dates[i-1]
        ncd = bond._flow_dates[i]
        if pcd < settlement_date and ncd > settlement_date:
            flow_time = (pcd - settlement_date) / gDaysInYear
            coupon_times.append(flow_time)
            coupon_flows.append(cpn)

    for flowDate in bond._flow_dates:
        if flowDate > settlement_date:
            flow_time = (flowDate - settlement_date) / gDaysInYear
            coupon_times.append(flow_time)
            coupon_flows.append(cpn)

    coupon_times = np.array(coupon_times)
    coupon_flows = np.array(coupon_flows)

    strikePrice = 105.0
    face = 100.0

    tmat = (maturity_date - settlement_date) / gDaysInYear
    texp = (expiry_date - settlement_date) / gDaysInYear
    times = np.linspace(0, tmat, 11)
    dates = settlement_date.addYears(times)
    dfs = np.exp(-0.05*times)

    testCases.header("LABEL", "VALUES")
    testCases.print("TIMES:", times)

    curve = DiscountCurve(settlement_date, dates, dfs)

    price = bond.clean_price_from_discount_curve(settlement_date, curve)
    testCases.print("Fixed Income Price:", price)

    sigma = 0.20

    # Test convergence
    num_stepsList = [5] #[100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    exerciseType = FinExerciseTypes.AMERICAN

    testCases.header("Values")
    treeVector = []
    for numTimeSteps in num_stepsList:
        model = FinModelRatesBDT(sigma, numTimeSteps)
        model.buildTree(tmat, times, dfs)
        v = model.bondOption(texp, strikePrice,
                             face, coupon_times, coupon_flows, exerciseType)

        testCases.print(v)
        treeVector.append(v['call'])

    if PLOT_GRAPHS:
        plt.plot(num_stepsList, treeVector)

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

    settlement_date = Date(1, 1, 2020)
    times = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    dates = settlement_date.addYears(times)
    rate = 0.06
    dfs = 1.0 / (1.0 + rate/2.0)**(2.0*times)
    curve = DiscountCurve(settlement_date, dates, dfs)

    coupon = 0.06
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    strikePrice = 100.0
    face = 100.0
    # Andersen paper
    numTimeSteps = 200

    testCases.header("ExerciseType", "Sigma", "NumSteps", "Texp", "Tmat", 
                     "V_Fixed", "V_pay", "V_rec")

    for exerciseType in [FinExerciseTypes.EUROPEAN,
                         FinExerciseTypes.BERMUDAN]:

        for maturityYears in [4.0, 5.0, 10.0, 20.0]:

            maturity_date = settlement_date.addYears(maturityYears)
            issue_date = Date(maturity_date._d, maturity_date._m, 2000)

            if maturityYears == 4.0 or maturityYears == 5.0:
                sigma = 0.2012
            elif maturityYears == 10.0:
                sigma = 0.1522
            elif maturityYears == 20.0:
                sigma = 0.1035

            for expiryYears in range(int(maturityYears/2)-1, int(maturityYears)):

                expiry_date = settlement_date.addYears(expiryYears)

                tmat = (maturity_date - settlement_date) / gDaysInYear
                texp = (expiry_date - settlement_date) / gDaysInYear

                bond = Bond(issue_date, maturity_date, coupon, freq_type, accrual_type)

                coupon_times = []
                coupon_flows = []
                cpn = bond._coupon/bond._frequency
                for flowDate in bond._flow_dates:
                    if flowDate > expiry_date:
                        flow_time = (flowDate - settlement_date) / gDaysInYear
                        coupon_times.append(flow_time)
                        coupon_flows.append(cpn)

                coupon_times = np.array(coupon_times)
                coupon_flows = np.array(coupon_flows)

                price = bond.clean_price_from_discount_curve(settlement_date, curve)

                model = FinModelRatesBDT(sigma, numTimeSteps)
                model.buildTree(tmat, times, dfs)

                v = model.bermudanSwaption(texp,
                                           tmat,
                                           strikePrice,
                                           face,
                                           coupon_times,
                                           coupon_flows,
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
