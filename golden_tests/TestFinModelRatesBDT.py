##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import sys
sys.path.append("..")

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
from financepy.utils.global_vars import g_days_in_year
from financepy.market.curves.discount_curve_zeros import DiscountCurveZeros
from financepy.models.bdt_tree import BDTTree
from financepy.utils.helpers import print_tree
from financepy.utils.global_types import FinExerciseTypes
from FinTestCases import FinTestCases, globalTestCaseMode


test_cases = FinTestCases(__file__, globalTestCaseMode)

PLOT_GRAPHS = False

###############################################################################


def testBlackModelCheck():

    # Checking Andersen paper using Black's model
    # Used to check swaption price below - we have Ts = 1 and Te = 4
    # Expect a price around 122 cents which is what I find.

    value_dt = Date(1, 1, 2020)
    libor_curve = DiscountCurveFlat(value_dt, 0.06,
                                    FrequencyTypes.SEMI_ANNUAL)

    settle_dt = Date(1, 1, 2020)
    exercise_dt = Date(1, 1, 2021)
    maturity_dt = Date(1, 1, 2024)

    fixed_cpn = 0.06
    fixed_freq_type = FrequencyTypes.SEMI_ANNUAL
    fixed_dc_type = DayCountTypes.THIRTY_E_360_ISDA
    notional = 100.0

    # Pricing a PAY
    swaptionType = SwapTypes.PAY
    swaption = IborSwaption(settle_dt,
                            exercise_dt,
                            maturity_dt,
                            swaptionType,
                            fixed_cpn,
                            fixed_freq_type,
                            fixed_dc_type,
                            notional)

    model = Black(0.20)
    v = swaption.value(value_dt, libor_curve, model)
    test_cases.header("LABEL", "VALUE")
    test_cases.print("BLACK'S MODEL PRICE:", v*100)

###############################################################################


def test_BDTExampleOne():
    # HULL BOOK NOTES
    # http://www-2.rotman.utoronto.ca/~hull/technicalnotes/TechnicalNote23.pdf

    value_dt = Date(1, 1, 2020)
    years = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    zero_dts = value_dt.add_years(years)
    zero_rates = [0.00, 0.10, 0.11, 0.12, 0.125, 0.13]

    test_cases.header("DATES")
    test_cases.print(zero_dts)

    test_cases.header("RATES")
    test_cases.print(zero_rates)

    curve = DiscountCurveZeros(value_dt,
                               zero_dts,
                               zero_rates,
                               FrequencyTypes.ANNUAL)

    yieldVol = 0.16

    num_time_steps = 5
    t_mat = years[-1]
    dfs = curve.df(zero_dts)

    test_cases.print("DFS")
    test_cases.print(dfs)

    years = np.array(years)
    dfs = np.array(dfs)

    model = BDTTree(yieldVol, num_time_steps)
    model.build_tree(t_mat, years, dfs)

###############################################################################


def test_BDTExampleTwo():
    # Valuation of a European option on a cpn bearing bond
    # This follows example in Fig 28.11 of John Hull's book (6th Edition)
    # but does not have the exact same dt so there are some differences

    test_cases.banner("===================== FIG 28.11 HULL BOOK =============")

    settle_dt = Date(1, 12, 2019)
    issue_dt = Date(1, 12, 2015)
    expiry_dt = settle_dt.add_tenor("18m")
    maturity_dt = settle_dt.add_tenor("10Y")
    cpn = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_dt, maturity_dt, cpn, freq_type, dc_type)

    cpn_times = []
    cpn_flows = []
    cpn = bond.cpn / bond.freq
    num_flows = len(bond.cpn_dts)

    for i in range(1, num_flows):
        pcd = bond.cpn_dts[i-1]
        ncd = bond.cpn_dts[i]
        if pcd < settle_dt and ncd > settle_dt:
            flow_time = (pcd - settle_dt) / g_days_in_year
            cpn_times.append(flow_time)
            cpn_flows.append(cpn)

    for flow_dt in bond.cpn_dts:
        if flow_dt > settle_dt:
            flow_time = (flow_dt - settle_dt) / g_days_in_year
            cpn_times.append(flow_time)
            cpn_flows.append(cpn)

    cpn_times = np.array(cpn_times)
    cpn_flows = np.array(cpn_flows)

    strike_price = 105.0
    face = 100.0

    t_mat = (maturity_dt - settle_dt) / g_days_in_year
    t_exp = (expiry_dt - settle_dt) / g_days_in_year
    times = np.linspace(0, t_mat, 11)
    dates = settle_dt.add_years(times)
    dfs = np.exp(-0.05*times)

    test_cases.header("LABEL", "VALUES")
    test_cases.print("TIMES:", times)

    curve = DiscountCurve(settle_dt, dates, dfs)

    price = bond.clean_price_from_discount_curve(settle_dt, curve)
    test_cases.print("Fixed Income Price:", price)

    sigma = 0.20

    # Test convergence
    num_steps_list = [5]  # [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    exercise_type = FinExerciseTypes.AMERICAN

    test_cases.header("Values")
    treeVector = []
    for num_time_steps in num_steps_list:
        model = BDTTree(sigma, num_time_steps)
        model.build_tree(t_mat, times, dfs)
        v = model.bond_option(t_exp, strike_price,
                              face, cpn_times, cpn_flows, exercise_type)

        test_cases.print(v)
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

    test_cases.banner("===================== ANDERSEN PAPER ==============")

    # This is a sanity check
    testBlackModelCheck()

    settle_dt = Date(1, 1, 2020)
    times = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    dates = settle_dt.add_years(times)
    rate = 0.06
    dfs = 1.0 / (1.0 + rate/2.0)**(2.0*times)
    curve = DiscountCurve(settle_dt, dates, dfs)

    cpn = 0.06
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    strike_price = 100.0
    face = 100.0
    # Andersen paper
    num_time_steps = 200

    test_cases.header("ExerciseType", "Sigma", "NumSteps", "Texp", "Tmat",
                     "V_Fixed", "V_pay", "V_rec")

    for exercise_type in [FinExerciseTypes.EUROPEAN,
                          FinExerciseTypes.BERMUDAN]:

        for years_to_maturity in [4.0, 5.0, 10.0, 20.0]:

            maturity_dt = settle_dt.add_years(years_to_maturity)
            issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)

            if years_to_maturity == 4.0 or years_to_maturity == 5.0:
                sigma = 0.2012
            elif years_to_maturity == 10.0:
                sigma = 0.1522
            elif years_to_maturity == 20.0:
                sigma = 0.1035

            for expiryYears in range(int(years_to_maturity/2)-1, int(years_to_maturity)):

                expiry_dt = settle_dt.add_years(expiryYears)

                t_mat = (maturity_dt - settle_dt) / g_days_in_year
                t_exp = (expiry_dt - settle_dt) / g_days_in_year

                bond = Bond(issue_dt, maturity_dt,
                            cpn, freq_type, dc_type)

                cpn_times = []
                cpn_flows = []
                cpn = bond.cpn / bond.freq
                for flow_dt in bond.cpn_dts:
                    if flow_dt > expiry_dt:
                        flow_time = (flow_dt - settle_dt) / g_days_in_year
                        cpn_times.append(flow_time)
                        cpn_flows.append(cpn)

                cpn_times = np.array(cpn_times)
                cpn_flows = np.array(cpn_flows)

                price = bond.clean_price_from_discount_curve(
                    settle_dt, curve)

                model = BDTTree(sigma, num_time_steps)
                model.build_tree(t_mat, times, dfs)

                v = model.bermudan_swaption(t_exp,
                                            t_mat,
                                            strike_price,
                                            face,
                                            cpn_times,
                                            cpn_flows,
                                            exercise_type)

                test_cases.print("%s" % exercise_type,
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

test_cases.compareTestCases()
