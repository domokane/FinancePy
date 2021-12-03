##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import matplotlib.pyplot as plt
import numpy as np
from financepy.utils.date import Date
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.products.bonds.bond import Bond
from financepy.products.rates.ibor_swaption import IborSwaption
from financepy.products.rates.ibor_swaption import SwapTypes
from financepy.models.black import Black
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.global_vars import gDaysInYear
from financepy.market.curves.discount_curve_zeros import DiscountCurveZeros
from financepy.models.bdt_tree import BDTTree
from financepy.utils.helpers import print_tree
from financepy.utils.global_types import FinExerciseTypes
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


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
    exercise_date = Date(1, 1, 2021)
    maturity_date = Date(1, 1, 2024)

    fixed_coupon = 0.06
    fixed_frequency_type = FrequencyTypes.SEMI_ANNUAL
    fixed_day_count_type = DayCountTypes.THIRTY_E_360_ISDA
    notional = 100.0

    # Pricing a PAY
    swaptionType = SwapTypes.PAY
    swaption = IborSwaption(settlement_date,
                            exercise_date,
                            maturity_date,
                            swaptionType,
                            fixed_coupon,
                            fixed_frequency_type,
                            fixed_day_count_type,
                            notional)

    model = Black(0.20)
    v = swaption.value(valuation_date, libor_curve, model)
    testCases.header("LABEL", "VALUE")
    testCases.print("BLACK'S MODEL PRICE:", v*100)

###############################################################################


def test_BDTExampleOne():
    # HULL BOOK NOTES
    # http://www-2.rotman.utoronto.ca/~hull/technicalnotes/TechnicalNote23.pdf

    valuation_date = Date(1, 1, 2020)
    years = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    zero_dates = valuation_date.add_years(years)
    zero_rates = [0.00, 0.10, 0.11, 0.12, 0.125, 0.13]

    testCases.header("DATES")
    testCases.print(zero_dates)

    testCases.header("RATES")
    testCases.print(zero_rates)

    curve = DiscountCurveZeros(valuation_date,
                               zero_dates,
                               zero_rates,
                               FrequencyTypes.ANNUAL)

    yieldVol = 0.16

    num_time_steps = 5
    tmat = years[-1]
    dfs = curve.df(zero_dates)

    testCases.print("DFS")
    testCases.print(dfs)

    years = np.array(years)
    dfs = np.array(dfs)

    model = BDTTree(yieldVol, num_time_steps)
    model.build_tree(tmat, years, dfs)

###############################################################################


def test_BDTExampleTwo():
    # Valuation of a European option on a coupon bearing bond
    # This follows example in Fig 28.11 of John Hull's book (6th Edition)
    # but does not have the exact same dt so there are some differences

    testCases.banner("===================== FIG 28.11 HULL BOOK =============")

    settlement_date = Date(1, 12, 2019)
    issue_date = Date(1, 12, 2015)
    expiry_date = settlement_date.add_tenor("18m")
    maturity_date = settlement_date.add_tenor("10Y")
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

    for flow_date in bond._flow_dates:
        if flow_date > settlement_date:
            flow_time = (flow_date - settlement_date) / gDaysInYear
            coupon_times.append(flow_time)
            coupon_flows.append(cpn)

    coupon_times = np.array(coupon_times)
    coupon_flows = np.array(coupon_flows)

    strike_price = 105.0
    face = 100.0

    tmat = (maturity_date - settlement_date) / gDaysInYear
    texp = (expiry_date - settlement_date) / gDaysInYear
    times = np.linspace(0, tmat, 11)
    dates = settlement_date.add_years(times)
    dfs = np.exp(-0.05*times)

    testCases.header("LABEL", "VALUES")
    testCases.print("TIMES:", times)

    curve = DiscountCurve(settlement_date, dates, dfs)

    price = bond.clean_price_from_discount_curve(settlement_date, curve)
    testCases.print("Fixed Income Price:", price)

    sigma = 0.20

    # Test convergence
    num_steps_list = [5]  # [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    exercise_type = FinExerciseTypes.AMERICAN

    testCases.header("Values")
    treeVector = []
    for num_time_steps in num_steps_list:
        model = BDTTree(sigma, num_time_steps)
        model.build_tree(tmat, times, dfs)
        v = model.bond_option(texp, strike_price,
                              face, coupon_times, coupon_flows, exercise_type)

        testCases.print(v)
        treeVector.append(v['call'])

    if PLOT_GRAPHS:
        plt.plot(num_steps_list, treeVector)

    # The value in Hull converges to 0.699 with 100 time steps while I get 0.70

    if 1 == 0:
        print("RT")
        print_tree(model._rt, 5)
        print("Q")
        print_tree(model._Q, 5)

###############################################################################


def test_BDTExampleThree():
    # Valuation of a swaption as in Leif Andersen's paper - see Table 1 on
    # SSRN-id155208.pdf

    testCases.banner("===================== ANDERSEN PAPER ==============")

    # This is a sanity check
    testBlackModelCheck()

    settlement_date = Date(1, 1, 2020)
    times = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    dates = settlement_date.add_years(times)
    rate = 0.06
    dfs = 1.0 / (1.0 + rate/2.0)**(2.0*times)
    curve = DiscountCurve(settlement_date, dates, dfs)

    coupon = 0.06
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    strike_price = 100.0
    face = 100.0
    # Andersen paper
    num_time_steps = 200

    testCases.header("ExerciseType", "Sigma", "NumSteps", "Texp", "Tmat",
                     "V_Fixed", "V_pay", "V_rec")

    for exercise_type in [FinExerciseTypes.EUROPEAN,
                          FinExerciseTypes.BERMUDAN]:

        for years_to_maturity in [4.0, 5.0, 10.0, 20.0]:

            maturity_date = settlement_date.add_years(years_to_maturity)
            issue_date = Date(maturity_date._d, maturity_date._m, 2000)

            if years_to_maturity == 4.0 or years_to_maturity == 5.0:
                sigma = 0.2012
            elif years_to_maturity == 10.0:
                sigma = 0.1522
            elif years_to_maturity == 20.0:
                sigma = 0.1035

            for expiryYears in range(int(years_to_maturity/2)-1, int(years_to_maturity)):

                expiry_date = settlement_date.add_years(expiryYears)

                tmat = (maturity_date - settlement_date) / gDaysInYear
                texp = (expiry_date - settlement_date) / gDaysInYear

                bond = Bond(issue_date, maturity_date,
                            coupon, freq_type, accrual_type)

                coupon_times = []
                coupon_flows = []
                cpn = bond._coupon/bond._frequency
                for flow_date in bond._flow_dates:
                    if flow_date > expiry_date:
                        flow_time = (flow_date - settlement_date) / gDaysInYear
                        coupon_times.append(flow_time)
                        coupon_flows.append(cpn)

                coupon_times = np.array(coupon_times)
                coupon_flows = np.array(coupon_flows)

                price = bond.clean_price_from_discount_curve(
                    settlement_date, curve)

                model = BDTTree(sigma, num_time_steps)
                model.build_tree(tmat, times, dfs)

                v = model.bermudan_swaption(texp,
                                            tmat,
                                            strike_price,
                                            face,
                                            coupon_times,
                                            coupon_flows,
                                            exercise_type)

                testCases.print("%s" % exercise_type,
                                "%9.5f" % sigma,
                                "%9.5f" % num_time_steps,
                                "%9.5f" % expiryYears,
                                "%9.5f" % years_to_maturity,
                                "%9.5f" % price,
                                "%9.2f" % (v['pay']*100.0),
                                "%9.2f" % (v['rec']*100.0))

###############################################################################
# This has broken and needs to be repaired!!!!


test_BDTExampleOne()
test_BDTExampleTwo()
test_BDTExampleThree()

testCases.compareTestCases()
