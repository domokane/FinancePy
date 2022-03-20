##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from financepy.utils.global_types import FinExerciseTypes
from financepy.utils.helpers import print_tree
from financepy.models.bdt_tree import BDTTree
from financepy.market.curves.discount_curve_zeros import DiscountCurveZeros
from financepy.utils.global_vars import gDaysInYear
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.models.black import Black
from financepy.products.rates.ibor_swaption import SwapTypes
from financepy.products.rates.ibor_swaption import IborSwaption
from financepy.products.bonds.bond import Bond
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.utils.date import Date
import numpy as np


def test_BDTExampleTwo():
    # Valuation of a European option on a coupon bearing bond
    # This follows example in Fig 28.11 of John Hull's book (6th Edition)
    # but does not have the exact same dt so there are some differences

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
    num_flows = len(bond._coupon_dates)

    for i in range(1, num_flows):
        pcd = bond._coupon_dates[i-1]
        ncd = bond._coupon_dates[i]
        if pcd < settlement_date and ncd > settlement_date:
            flow_time = (pcd - settlement_date) / gDaysInYear
            coupon_times.append(flow_time)
            coupon_flows.append(cpn)

    for flow_date in bond._coupon_dates:
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

    curve = DiscountCurve(settlement_date, dates, dfs)

    price = bond.clean_price_from_discount_curve(settlement_date, curve)
    assert round(price, 4) == 99.5420

    sigma = 0.20

    # Test convergence
    num_time_steps = 5
    exercise_type = FinExerciseTypes.AMERICAN

    model = BDTTree(sigma, num_time_steps)
    model.build_tree(tmat, times, dfs)
    v = model.bond_option(texp, strike_price,
                          face, coupon_times, coupon_flows, exercise_type)

    assert round(v['call'], 4) == 0.5043
    assert round(v['put'], 4) == 8.2242


def test_BDTExampleThree():
    # Valuation of a swaption as in Leif Andersen's paper - see Table 1 on
    # SSRN-id155208.pdf

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

    exercise_type = FinExerciseTypes.EUROPEAN
    years_to_maturity = 4.0
    expiryYears = 2.0

    maturity_date = settlement_date.add_years(years_to_maturity)
    issue_date = Date(maturity_date._d, maturity_date._m, 2000)

    sigma = 0.2012

    expiry_date = settlement_date.add_years(expiryYears)

    tmat = (maturity_date - settlement_date) / gDaysInYear
    texp = (expiry_date - settlement_date) / gDaysInYear

    bond = Bond(issue_date, maturity_date,
                coupon, freq_type, accrual_type)

    coupon_times = []
    coupon_flows = []
    cpn = bond._coupon/bond._frequency
    for flow_date in bond._coupon_dates:
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

    assert round(price, 5) == 100.01832
    assert round(v['pay']*100, 2) == 0.00
    assert round(v['rec']*100, 2) == 8883.21

    exercise_type = FinExerciseTypes.BERMUDAN
    years_to_maturity = 10.0
    expiryYears = 5.0

    maturity_date = settlement_date.add_years(years_to_maturity)
    issue_date = Date(maturity_date._d, maturity_date._m, 2000)

    sigma = 0.1522

    expiry_date = settlement_date.add_years(expiryYears)

    tmat = (maturity_date - settlement_date) / gDaysInYear
    texp = (expiry_date - settlement_date) / gDaysInYear

    bond = Bond(issue_date, maturity_date,
                coupon, freq_type, accrual_type)

    coupon_times = []
    coupon_flows = []
    cpn = bond._coupon/bond._frequency
    for flow_date in bond._coupon_dates:
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

    assert round(price, 5) == 100.08625
    assert round(v['pay']*100, 2) == 263.28
    assert round(v['rec']*100, 2) == 7437.00
