# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
import pytest

from financepy.utils.error import FinError
from financepy.utils.global_types import ExerciseTypes
from financepy.models.bdt_tree import BDTTree
from financepy.utils.global_vars import G_DAYS_IN_YEAR
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.bonds.bond import Bond
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.utils.date import Date

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
            flow_time = (pcd - settle_dt) / G_DAYS_IN_YEAR
            cpn_times.append(flow_time)
            cpn_flows.append(cpn)

    for flow_dt in bond.cpn_dts:
        if flow_dt > settle_dt:
            flow_time = (flow_dt - settle_dt) / G_DAYS_IN_YEAR
            cpn_times.append(flow_time)
            cpn_flows.append(cpn)

    cpn_times = np.array(cpn_times)
    cpn_flows = np.array(cpn_flows)

    strike_price = 105.0
    face = 100.0

    t_mat = (maturity_dt - settle_dt) / G_DAYS_IN_YEAR
    t_exp = (expiry_dt - settle_dt) / G_DAYS_IN_YEAR
    times = np.linspace(0, t_mat, 11)
    dates = settle_dt.add_years(times)
    dfs = np.exp(-0.05 * times)

    curve = DiscountCurve(settle_dt, dates, dfs)

    price = bond.clean_price_from_discount_curve(settle_dt, curve)
    assert round(price, 4) == 99.4906

    sigma = 0.20

    # Test convergence
    num_time_steps = 5
    exercise_type = ExerciseTypes.AMERICAN

    model = BDTTree(sigma, num_time_steps)
    model.build_tree(t_mat, times, dfs)
    v = model.bond_option(
        t_exp, strike_price, face, cpn_times, cpn_flows, exercise_type
    )

    assert round(v[0], 4) == 0.6121
    assert round(v[1], 4) == 7.9605


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

    exercise_type = ExerciseTypes.EUROPEAN
    years_to_maturity = 4.0
    expiry_years = 2.0

    maturity_dt = settle_dt.add_years(years_to_maturity)
    issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)

    sigma = 0.2012

    expiry_dt = settle_dt.add_years(expiry_years)

    t_mat = (maturity_dt - settle_dt) / G_DAYS_IN_YEAR
    t_exp = (expiry_dt - settle_dt) / G_DAYS_IN_YEAR

    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    cpn_times = []
    cpn_flows = []

    cpn = bond.cpn / bond.freq

    for flow_dt in bond.cpn_dts:
        if flow_dt > expiry_dt:
            flow_time = (flow_dt - settle_dt) / G_DAYS_IN_YEAR
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

    assert round(price, 5) == 100.0003
    assert round(v[0] * 100, 2) == 0.00
    assert round(v[1] * 100, 2) == 8883.55

    exercise_type = ExerciseTypes.BERMUDAN
    years_to_maturity = 10.0
    expiry_years = 5.0

    maturity_dt = settle_dt.add_years(years_to_maturity)
    issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)

    sigma = 0.1522

    expiry_dt = settle_dt.add_years(expiry_years)

    t_mat = (maturity_dt - settle_dt) / G_DAYS_IN_YEAR
    t_exp = (expiry_dt - settle_dt) / G_DAYS_IN_YEAR

    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    cpn_times = []
    cpn_flows = []

    cpn = bond.cpn / bond.freq

    for flow_dt in bond.cpn_dts:
        if flow_dt > expiry_dt:
            flow_time = (flow_dt - settle_dt) / G_DAYS_IN_YEAR
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

    assert round(price, 5) == 100.04429
    assert round(v[0] * 100, 2) == 263.27
    assert round(v[1] * 100, 2) == 7437.75


########################################################################################


def test_shifted_bdt_negative_rates_calibration_and_spacing():
    """Shifted BDT tree calibrates to negative-rate discount curves and
    preserves lognormal spacing on the shifted short rate (r + shift)."""
    t_mat = 5.0
    num_time_steps = 20
    sigma = 0.20
    shift = 0.04
    neg_rate = -0.015

    times = np.linspace(0.0, t_mat, 11)
    dfs = np.exp(-neg_rate * times)

    model = BDTTree(sigma, num_time_steps, shift=shift)
    model.build_tree(t_mat, times, dfs)

    assert model.shift == pytest.approx(shift)
    assert model.rt[0, 0] == pytest.approx(neg_rate, abs=1e-12)

    dt = model._dt
    expected_step = 2.0 * sigma * np.sqrt(dt)

    for m in range(0, num_time_steps + 2):
        # Arrow-Debreu state prices must recover the discount factor at each step
        q_sum = np.sum(model.qq[m, : m + 1])
        expected_df = np.exp(-neg_rate * model._tree_times[m])
        assert q_sum == pytest.approx(expected_df, abs=1e-9)

    for m in range(1, num_time_steps + 1):
        shifted_rates = model.rt[m, : m + 1] + shift
        assert np.all(shifted_rates > 0.0)
        log_ratios = np.log(shifted_rates[1:] / shifted_rates[:-1])
        np.testing.assert_allclose(log_ratios, expected_step, rtol=1e-12, atol=1e-12)

    # Put-call parity for European bond option under negative interest rates
    t_exp = 2.0
    cpn_times = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    cpn_flows = np.array([0.02, 0.02, 0.02, 0.02, 0.02])
    call_v, put_v = model.bond_option(
        t_exp, 100.0, 100.0, cpn_times, cpn_flows, ExerciseTypes.EUROPEAN
    )
    fwd_bond_pv = (
        0.02 * np.exp(-neg_rate * 3.0)
        + 0.02 * np.exp(-neg_rate * 4.0)
        + 1.02 * np.exp(-neg_rate * 5.0)
    ) * 100.0
    strike_pv = 100.0 * np.exp(-neg_rate * t_exp)
    assert (call_v - put_v) == pytest.approx(fwd_bond_pv - strike_pv, abs=1e-8)


def test_shifted_bdt_zero_vol_shift_invariance_and_pricing():
    """At zero volatility, shifted BDT short rates and option prices on a
    negative-rate curve are invariant to the choice of valid shift."""
    t_mat = 4.0
    t_exp = 2.0
    num_time_steps = 40
    neg_rate = -0.01

    times = np.linspace(0.0, t_mat, 9)
    dfs = np.exp(-neg_rate * times)

    cpn_times = np.array([1.0, 2.0, 3.0, 4.0])
    cpn_flows = np.array([0.01, 0.01, 0.01, 0.01])

    model_s1 = BDTTree(0.0, num_time_steps, shift=0.03)
    model_s2 = BDTTree(0.0, num_time_steps, shift=0.08)

    model_s1.build_tree(t_mat, times, dfs)
    model_s2.build_tree(t_mat, times, dfs)

    for m in range(0, num_time_steps + 1):
        np.testing.assert_allclose(
            model_s1.rt[m, : m + 1],
            model_s2.rt[m, : m + 1],
            atol=1e-10,
        )

    v_bond_s1 = model_s1.bond_option(
        t_exp, 100.0, 100.0, cpn_times, cpn_flows, ExerciseTypes.EUROPEAN
    )
    v_bond_s2 = model_s2.bond_option(
        t_exp, 100.0, 100.0, cpn_times, cpn_flows, ExerciseTypes.EUROPEAN
    )
    np.testing.assert_allclose(v_bond_s1, v_bond_s2, atol=1e-8)

    v_swap_s1 = model_s1.bermudan_swaption(
        t_exp, t_mat, 100.0, 100.0, cpn_times[1:], cpn_flows[1:], ExerciseTypes.BERMUDAN
    )
    v_swap_s2 = model_s2.bermudan_swaption(
        t_exp, t_mat, 100.0, 100.0, cpn_times[1:], cpn_flows[1:], ExerciseTypes.BERMUDAN
    )
    np.testing.assert_allclose(v_swap_s1, v_swap_s2, atol=1e-8)


def test_shifted_bdt_rejects_rates_below_negative_shift():
    """Building a BDT tree must raise FinError if r + shift <= 0."""
    t_mat = 3.0
    times = np.array([0.0, 1.0, 2.0, 3.0])
    dfs_neg = np.exp(0.02 * times)  # r = -2%

    # Unshifted BDT cannot fit negative initial short rate
    model_unshifted = BDTTree(0.15, 10)
    with pytest.raises(FinError, match="Short rate plus shift must be positive"):
        model_unshifted.build_tree(t_mat, times, dfs_neg)

    # Insufficient shift (shift = 1% for r = -2%) must also raise FinError
    model_small_shift = BDTTree(0.15, 10, shift=0.01)
    with pytest.raises(FinError, match="Short rate plus shift must be positive"):
        model_small_shift.build_tree(t_mat, times, dfs_neg)

    # Curve starting positive (+1%) then dropping to -3% with shift = 2%
    dfs_later_neg = np.array([1.0, np.exp(-0.01), np.exp(-0.01 + 0.03), np.exp(-0.01 + 0.06)])
    model_mid_fail = BDTTree(0.15, 6, shift=0.02)
    with pytest.raises(FinError, match="Short rate plus shift must be positive"):
        model_mid_fail.build_tree(t_mat, times, dfs_later_neg)
