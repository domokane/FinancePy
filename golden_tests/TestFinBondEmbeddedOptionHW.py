###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import matplotlib.pyplot as plt
import time

import sys

sys.path.append("..")

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.bonds.bond_callable import BondEmbeddedOption
from financepy.products.bonds.bond import Bond
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_swap import IborSwap
from financepy.models.hw_tree import HWTree
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date
from financepy.utils.global_types import SwapTypes

test_cases = FinTestCases(__file__, globalTestCaseMode)

plotGraphs = False

###############################################################################


def test_BondEmbeddedOptionMATLAB():

    # https://fr.mathworks.com/help/fininst/optembndbyhw.html
    # I FIND THAT THE PRICE CONVERGES TO 102.88 WHICH IS CLOSE TO 102.9127
    # FOUND BY MATLAB ALTHOUGH THEY DO NOT EXAMINE THE ASYMPTOTIC PRICE
    # WHICH MIGHT BE A BETTER MATCH

    settle_dt = Date(1, 1, 2007)
    value_dt = settle_dt

    ###########################################################################

    dc_type = DayCountTypes.THIRTY_E_360
    fixed_freq = FrequencyTypes.ANNUAL
    fixed_leg_type = SwapTypes.PAY
    swap1 = IborSwap(
        settle_dt, "1Y", fixed_leg_type, 0.0350, fixed_freq, dc_type
    )
    swap2 = IborSwap(
        settle_dt, "2Y", fixed_leg_type, 0.0400, fixed_freq, dc_type
    )
    swap3 = IborSwap(
        settle_dt, "3Y", fixed_leg_type, 0.0450, fixed_freq, dc_type
    )
    swaps = [swap1, swap2, swap3]
    discount_curve = IborSingleCurve(value_dt, [], [], swaps)

    ###########################################################################

    issue_dt = Date(1, 1, 2004)
    maturity_dt = Date(1, 1, 2010)

    coupon = 0.0525
    freq_type = FrequencyTypes.ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    call_dts = []
    call_prices = []
    put_dts = []
    put_prices = []

    put_dt = Date(1, 1, 2008)
    for _ in range(0, 24):
        put_dts.append(put_dt)
        put_prices.append(100)
        put_dt = put_dt.add_months(1)

    test_cases.header("BOND PRICE", "PRICE")
    v = bond.clean_price_from_discount_curve(settle_dt, discount_curve)
    test_cases.print("Bond Pure Price:", v)

    sigma = 0.01  # basis point volatility
    a = 0.1

    puttable_bond = BondEmbeddedOption(
        issue_dt,
        maturity_dt,
        coupon,
        freq_type,
        dc_type,
        call_dts,
        call_prices,
        put_dts,
        put_prices,
    )

    test_cases.header("TIME", "Numtime_steps", "BondWithOption", "BondPure")

    time_steps = range(50, 1000, 50)
    values = []
    for num_time_steps in time_steps:
        model = HWTree(sigma, a, num_time_steps)
        start = time.time()
        v = puttable_bond.value(settle_dt, discount_curve, model)
        end = time.time()
        period = end - start
        test_cases.print(
            period, num_time_steps, v["bondwithoption"], v["bondpure"]
        )
        values.append(v["bondwithoption"])

    if plotGraphs:
        plt.figure()
        plt.plot(time_steps, values)


###############################################################################


def test_BondEmbeddedOptionQUANTLIB():

    # Based on example at the nice blog on Quantlib at
    # http://gouthamanbalaraman.com/blog/callable-bond-quantlib-python.html
    # I get a price of 68.97 for 1000 time steps which is higher than the
    # 68.38 found in blog article. But this is for 40 grid points.
    # Note also that a basis point vol of 0.120 is 12% which is VERY HIGH!

    value_dt = Date(16, 8, 2016)
    settle_dt = value_dt.add_weekdays(3)

    ###########################################################################

    discount_curve = DiscountCurveFlat(
        value_dt, 0.035, FrequencyTypes.SEMI_ANNUAL
    )

    ###########################################################################

    issue_dt = Date(15, 9, 2010)
    maturity_dt = Date(15, 9, 2022)
    coupon = 0.025
    freq_type = FrequencyTypes.QUARTERLY
    dc_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    ###########################################################################
    # Set up the call and put times and prices
    ###########################################################################

    nextCallDate = Date(15, 9, 2016)
    call_dts = [nextCallDate]
    call_prices = [100.0]

    for _ in range(1, 24):
        nextCallDate = nextCallDate.add_months(3)
        call_dts.append(nextCallDate)
        call_prices.append(100.0)

    put_dts = []
    put_prices = []

    # the value used in blog of 12% bp vol is unrealistic
    sigma = 0.12  # basis point volatility
    a = 0.03

    puttable_bond = BondEmbeddedOption(
        issue_dt,
        maturity_dt,
        coupon,
        freq_type,
        dc_type,
        call_dts,
        call_prices,
        put_dts,
        put_prices,
    )

    test_cases.header("BOND PRICE", "PRICE")
    v = bond.clean_price_from_discount_curve(settle_dt, discount_curve)
    test_cases.print("Bond Pure Price:", v)

    test_cases.header("TIME", "Numtime_steps", "BondWithOption", "BondPure")
    time_steps = range(100, 200, 50)
    values = []
    for num_time_steps in time_steps:
        model = HWTree(sigma, a, num_time_steps)
        start = time.time()
        v = puttable_bond.value(settle_dt, discount_curve, model)
        end = time.time()
        period = end - start
        test_cases.print(
            period, num_time_steps, v["bondwithoption"], v["bondpure"]
        )
        values.append(v["bondwithoption"])

    if plotGraphs:
        plt.figure()
        plt.title("Puttable Bond Price Convergence")
        plt.plot(time_steps, values)


###############################################################################


test_BondEmbeddedOptionMATLAB()
test_BondEmbeddedOptionQUANTLIB()
test_cases.compareTestCases()
