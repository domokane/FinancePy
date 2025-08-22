# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.utils.global_types import FinExerciseTypes
from financepy.utils.helpers import print_tree
from financepy.models.bdt_tree import BDTTree
from financepy.market.curves.discount_curve_zeros import DiscountCurveZeros
from financepy.utils.global_vars import G_DAYS_IN_YEARS
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

########################################################################################


def test_bdt_example_two():

    # Valuation of a European option on a coupon bearing bond
    # This follows example in Fig 28.11 of John Hull's book (6th Edition)
    # but does not have the exact same dt so there are some differences

    settle_dt = Date(1, 12, 2019)
    issue_dt = Date(1, 12, 2015)
    expiry_dt = settle_dt.add_tenor("18m")
    maturity_dt = settle_dt.add_tenor("10Y")
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    cpn_times = []
    cpn_flows = []
    cpn = bond.cpn / bond.freq
    num_flows = len(bond.cpn_dts)

    for i in range(1, num_flows):
        pcd = bond.cpn_dts[i - 1]
        ncd = bond.cpn_dts[i]
        if pcd < settle_dt and ncd > settle_dt:
            flow_time = (pcd - settle_dt) / G_DAYS_IN_YEARS
            cpn_times.append(flow_time)
            cpn_flows.append(cpn)

    for flow_dt in bond.cpn_dts:
        if flow_dt > settle_dt:
            flow_time = (flow_dt - settle_dt) / G_DAYS_IN_YEARS
            cpn_times.append(flow_time)
            cpn_flows.append(cpn)

    cpn_times = np.array(cpn_times)
    cpn_flows = np.array(cpn_flows)

    strike_price = 105.0
    face = 100.0

    t_mat = (maturity_dt - settle_dt) / G_DAYS_IN_YEARS
    t_exp = (expiry_dt - settle_dt) / G_DAYS_IN_YEARS
    times = np.linspace(0, t_mat, 11)
    dates = settle_dt.add_years(times)
    dfs = np.exp(-0.05 * times)

    curve = DiscountCurve(settle_dt, dates, dfs)

    price = bond.clean_price_from_discount_curve(settle_dt, curve)
    assert round(price, 4) == 99.5420

    sigma = 0.20

    # Test convergence
    num_time_steps = 5
    exercise_type = FinExerciseTypes.AMERICAN

    model = BDTTree(sigma, num_time_steps)
    model.build_tree(t_mat, times, dfs)
    v = model.bond_option(
        t_exp, strike_price, face, cpn_times, cpn_flows, exercise_type
    )

    assert round(v["call"], 4) == 0.5043
    assert round(v["put"], 4) == 8.2242

########################################################################################


def test_bdt_example_three():

    # Valuation of a swaption as in Leif Andersen's paper - see Table 1 on
    # SSRN-id155208.pdf

    settle_dt = Date(1, 1, 2020)
    times = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    dates = settle_dt.add_years(times)
    rate = 0.06
    dfs = 1.0 / (1.0 + rate / 2.0) ** (2.0 * times)
    curve = DiscountCurve(settle_dt, dates, dfs)

    coupon = 0.06
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    strike_price = 100.0
    face = 100.0
    # Andersen paper
    num_time_steps = 200

    exercise_type = FinExerciseTypes.EUROPEAN
    years_to_maturity = 4.0
    expiry_years = 2.0

    maturity_dt = settle_dt.add_years(years_to_maturity)
    issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)

    sigma = 0.2012

    expiry_dt = settle_dt.add_years(expiry_years)

    t_mat = (maturity_dt - settle_dt) / G_DAYS_IN_YEARS
    t_exp = (expiry_dt - settle_dt) / G_DAYS_IN_YEARS

    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    cpn_times = []
    cpn_flows = []

    cpn = bond.cpn / bond.freq

    for flow_dt in bond.cpn_dts:
        if flow_dt > expiry_dt:
            flow_time = (flow_dt - settle_dt) / G_DAYS_IN_YEARS
            cpn_times.append(flow_time)
            cpn_flows.append(cpn)

    cpn_times = np.array(cpn_times)
    cpn_flows = np.array(cpn_flows)

    price = bond.clean_price_from_discount_curve(settle_dt, curve)

    model = BDTTree(sigma, num_time_steps)
    model.build_tree(t_mat, times, dfs)

    v = model.bermudan_swaption(
        t_exp, t_mat, strike_price, face, cpn_times, cpn_flows, exercise_type
    )

    assert round(price, 5) == 100.01832
    assert round(v["pay"] * 100, 2) == 0.00
    assert round(v["rec"] * 100, 2) == 8883.21

    exercise_type = FinExerciseTypes.BERMUDAN
    years_to_maturity = 10.0
    expiry_years = 5.0

    maturity_dt = settle_dt.add_years(years_to_maturity)
    issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)

    sigma = 0.1522

    expiry_dt = settle_dt.add_years(expiry_years)

    t_mat = (maturity_dt - settle_dt) / G_DAYS_IN_YEARS
    t_exp = (expiry_dt - settle_dt) / G_DAYS_IN_YEARS

    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    cpn_times = []
    cpn_flows = []

    cpn = bond.cpn / bond.freq

    for flow_dt in bond.cpn_dts:
        if flow_dt > expiry_dt:
            flow_time = (flow_dt - settle_dt) / G_DAYS_IN_YEARS
            cpn_times.append(flow_time)
            cpn_flows.append(cpn)

    cpn_times = np.array(cpn_times)
    cpn_flows = np.array(cpn_flows)

    price = bond.clean_price_from_discount_curve(settle_dt, curve)

    model = BDTTree(sigma, num_time_steps)
    model.build_tree(t_mat, times, dfs)

    v = model.bermudan_swaption(
        t_exp, t_mat, strike_price, face, cpn_times, cpn_flows, exercise_type
    )

    assert round(price, 5) == 100.08625
    assert round(v["pay"] * 100, 2) == 263.28
    assert round(v["rec"] * 100, 2) == 7437.00
