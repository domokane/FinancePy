###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

import time

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.models.bk_tree import BKTree
from financepy.utils.global_types import SwapTypes
from financepy.products.bonds.bond_callable import BondEmbeddedOption
from financepy.products.bonds.bond import Bond
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_swap import IborSwap
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date
import matplotlib.pyplot as plt


testCases = FinTestCases(__file__, globalTestCaseMode)

plotGraphs = False

###############################################################################


def test_BondEmbeddedOptionMATLAB():
    # https://fr.mathworks.com/help/fininst/optembndbybk.html
    # I FIND THAT THE PRICE CONVERGES TO 102.365 WHICH IS CLOSE TO 102.382
    # FOUND BY MATLAB ALTHOUGH THEY DO NOT EXAMINE THE ASYMPTOTIC PRICE
    # WHICH MIGHT BE A BETTER MATCH - ALSO THEY DO NOT USE A REALISTIC VOL

    valuation_date = Date(1, 1, 2007)
    settlement_date = valuation_date

    ###########################################################################

    fixed_leg_type = SwapTypes.PAY
    dcType = DayCountTypes.THIRTY_E_360
    fixedFreq = FrequencyTypes.ANNUAL
    swap1 = IborSwap(settlement_date, "1Y", fixed_leg_type,
                     0.0350, fixedFreq, dcType)
    swap2 = IborSwap(settlement_date, "2Y", fixed_leg_type,
                     0.0400, fixedFreq, dcType)
    swap3 = IborSwap(settlement_date, "3Y", fixed_leg_type,
                     0.0450, fixedFreq, dcType)
    swaps = [swap1, swap2, swap3]
    discount_curve = IborSingleCurve(valuation_date, [], [], swaps)

    ###########################################################################

    issue_date = Date(1, 1, 2005)
    maturity_date = Date(1, 1, 2010)
    coupon = 0.0525
    freq_type = FrequencyTypes.ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_date, maturity_date, coupon, freq_type, accrual_type)

    call_dates = []
    call_prices = []
    put_dates = []
    put_prices = []

    putDate = Date(1, 1, 2008)
    for _ in range(0, 24):
        put_dates.append(putDate)
        put_prices.append(100)
        putDate = putDate.add_months(1)

    testCases.header("BOND PRICE", "PRICE")
    v = bond.clean_price_from_discount_curve(settlement_date, discount_curve)
    testCases.print("Bond Pure Price:", v)

    sigma = 0.01  # This volatility is very small for a BK process
    a = 0.1

    puttableBond = BondEmbeddedOption(issue_date, maturity_date, coupon,
                                      freq_type, accrual_type,
                                      call_dates, call_prices,
                                      put_dates, put_prices)

    testCases.header("TIME", "NumTimeSteps", "BondWithOption", "BondPure")

    timeSteps = range(100, 200, 10)  # 1000, 10)
    values = []
    for num_time_steps in timeSteps:
        model = BKTree(sigma, a, num_time_steps)
        start = time.time()
        v = puttableBond.value(settlement_date, discount_curve, model)
        end = time.time()
        period = end - start
        testCases.print(period, num_time_steps, v['bondwithoption'],
                        v['bondpure'])

        values.append(v['bondwithoption'])

    if plotGraphs:
        plt.figure()
        plt.plot(timeSteps, values)

###############################################################################


def test_BondEmbeddedOptionQUANTLIB():

    # Based on example at the nice blog on Quantlib at
    # http://gouthamanbalaraman.com/blog/callable-bond-quantlib-python.html
    # I get a price of 68.97 for 1000 time steps which is higher than the
    # 68.38 found in blog article. But this is for 40 grid points.
    # Note also that a basis point vol of 0.120 is 12% which is VERY HIGH!

    valuation_date = Date(16, 8, 2016)
    settlement_date = valuation_date.add_weekdays(3)

    ###########################################################################

    discount_curve = DiscountCurveFlat(valuation_date, 0.035,
                                       FrequencyTypes.SEMI_ANNUAL)

    ###########################################################################

    issue_date = Date(15, 9, 2010)
    maturity_date = Date(15, 9, 2022)
    coupon = 0.025
    freq_type = FrequencyTypes.QUARTERLY
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_date, maturity_date, coupon, freq_type, accrual_type)

    ###########################################################################
    # Set up the call and put times and prices
    ###########################################################################

    nextCallDate = Date(15, 9, 2016)
    call_dates = [nextCallDate]
    call_prices = [100.0]

    for _ in range(1, 24):
        nextCallDate = nextCallDate.add_months(3)
        call_dates.append(nextCallDate)
        call_prices.append(100.0)

    put_dates = []
    put_prices = []

    # the value used in blog of 12% bp vol is unrealistic
    sigma = 0.12/0.035  # basis point volatility
    a = 0.03

    puttableBond = BondEmbeddedOption(issue_date, maturity_date, coupon,
                                      freq_type, accrual_type,
                                      call_dates, call_prices,
                                      put_dates, put_prices)

    testCases.header("BOND PRICE", "PRICE")
    v = bond.clean_price_from_discount_curve(settlement_date, discount_curve)
    testCases.print("Bond Pure Price:", v)

    testCases.header("TIME", "NumTimeSteps", "BondWithOption", "BondPure")
    timeSteps = range(100, 200, 20)  # 1000, 10)
    values = []
    for num_time_steps in timeSteps:
        model = BKTree(sigma, a, num_time_steps)
        start = time.time()
        v = puttableBond.value(settlement_date, discount_curve, model)
        end = time.time()
        period = end - start
        testCases.print(period, num_time_steps, v['bondwithoption'],
                        v['bondpure'])
        values.append(v['bondwithoption'])

    if plotGraphs:
        plt.figure()
        plt.title("Puttable Bond Price Convergence")
        plt.plot(timeSteps, values)

###############################################################################


test_BondEmbeddedOptionMATLAB()
test_BondEmbeddedOptionQUANTLIB()
testCases.compareTestCases()
