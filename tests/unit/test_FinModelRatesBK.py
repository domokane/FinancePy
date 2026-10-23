# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
import pytest

from financepy.utils.error import FinError

from financepy.utils.global_types import ExerciseTypes
from financepy.models.bk_tree import BKTree
from financepy.utils.global_vars import G_DAYS_IN_YEAR
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.bonds.bond import Bond
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.utils.date import Date

########################################################################################


def test_bk_example_one():

    times = [0.0, 0.5000, 1.00000, 1.50000, 2.00000, 2.500000, 3.00000]
    zeros = [0.03, 0.0343, 0.03824, 0.04183, 0.04512, 0.048512, 0.05086]
    times = np.array(times)
    zeros = np.array(zeros)
    dfs = np.exp(-zeros * times)

    start_dt = Date(1, 12, 2019)
    end_dt = Date(1, 6, 2021)
    sigma = 0.25
    a = 0.22
    num_time_steps = 3
    t_mat = (end_dt - start_dt) / G_DAYS_IN_YEAR
    model = BKTree(sigma, a, num_time_steps)
    model.build_tree(t_mat, times, dfs)

    assert [round(x, 4) for x in model.qq[2]] == [
        0.0190,
        0.2126,
        0.5009,
        0.2112,
        0.0187,
    ]

    assert [round(x, 4) for x in model.rt[2]] == [
        0.0259,
        0.0351,
        0.0477,
        0.0648,
        0.0881,
    ]

    assert [round(x, 4) for x in model.pu] == [
        0.0808,
        0.2278,
        0.1667,
        0.1177,
        0.8606,
    ]

    assert [round(x, 4) for x in model.pm] == [
        0.0586,
        0.6545,
        0.6667,
        0.6545,
        0.0586,
    ]

    assert [round(x, 4) for x in model.pd] == [
        0.8606,
        0.1177,
        0.1667,
        0.2278,
        0.0808,
    ]


########################################################################################


def test_bk_example_two():

    # Valuation of a European option on a coupon bearing bond
    # This follows example in Fig 28.11 of John Hull's book but does not
    # have the exact same dt so there are some differences

    settle_dt = Date(1, 12, 2019)
    issue_dt = Date(1, 12, 2018)
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
    a = 0.05
    num_time_steps = 26

    model = BKTree(sigma, a, num_time_steps)
    model.build_tree(t_mat, times, dfs)
    exercise_type = ExerciseTypes.AMERICAN
    v = model.bond_option(
        t_exp, strike_price, face, cpn_times, cpn_flows, exercise_type
    )

    # Test convergence
    num_time_steps = 200
    exercise_type = ExerciseTypes.AMERICAN

    tree_vector = []
    model = BKTree(sigma, a, num_time_steps)
    model.build_tree(t_mat, times, dfs)
    v = model.bond_option(
        t_exp, strike_price, face, cpn_times, cpn_flows, exercise_type
    )
    tree_vector.append(v)

    assert round(v[0], 4) == 0.6998
    assert round(v[1], 4) == 7.9605


########################################################################################

test_bk_example_one()

########################################################################################


def test_shifted_bk_negative_rates_calibration_and_spacing():
    """Shifted BK tree calibrates to negative-rate discount curves and
    preserves lognormal spacing on the shifted short rate (r + shift)."""
    t_mat = 5.0
    num_time_steps = 20
    sigma = 0.20
    a = 0.05
    shift = 0.04
    neg_rate = -0.015

    times = np.linspace(0.0, t_mat, 11)
    dfs = np.exp(-neg_rate * times)

    model = BKTree(sigma, a, num_time_steps, shift=shift)
    model.build_tree(t_mat, times, dfs)

    assert model.shift == pytest.approx(shift)
    j_max = len(model.pu) // 2
    assert model.rt[0, j_max] == pytest.approx(neg_rate, abs=1e-12)

    dt = model.dt
    expected_step = sigma * np.sqrt(3.0 * dt)

    for m in range(0, num_time_steps + 2):
        # Arrow-Debreu state prices must recover the discount factor at each step
        q_sum = np.sum(model.qq[m])
        expected_df = np.exp(-neg_rate * model.tree_times[m])
        assert q_sum == pytest.approx(expected_df, abs=1e-9)

    for m in range(1, num_time_steps + 1):
        nm = min(m, j_max)
        shifted_rates = model.rt[m, j_max - nm : j_max + nm + 1] + shift
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


def test_shifted_bk_zero_vol_shift_invariance_and_pricing():
    """At zero volatility, shifted BK short rates and option prices on a
    negative-rate curve are invariant to the choice of valid shift."""
    t_mat = 4.0
    t_exp = 2.0
    num_time_steps = 40
    neg_rate = -0.01

    times = np.linspace(0.0, t_mat, 9)
    dfs = np.exp(-neg_rate * times)

    cpn_times = np.array([1.0, 2.0, 3.0, 4.0])
    cpn_flows = np.array([0.01, 0.01, 0.01, 0.01])

    # Zero vol with mean reversion
    model_s1 = BKTree(0.0, 0.05, num_time_steps, shift=0.03)
    model_s2 = BKTree(0.0, 0.05, num_time_steps, shift=0.08)

    model_s1.build_tree(t_mat, times, dfs)
    model_s2.build_tree(t_mat, times, dfs)

    j_max_s1 = len(model_s1.pu) // 2
    j_max_s2 = len(model_s2.pu) // 2
    assert j_max_s1 == j_max_s2

    for m in range(0, num_time_steps + 1):
        nm = min(m, j_max_s1)
        np.testing.assert_allclose(
            model_s1.rt[m, j_max_s1 - nm : j_max_s1 + nm + 1],
            model_s2.rt[m, j_max_s2 - nm : j_max_s2 + nm + 1],
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


def test_shifted_bk_rejects_rates_below_negative_shift():
    """Building a BK tree must raise FinError if r + shift <= 0."""
    t_mat = 3.0
    times = np.array([0.0, 1.0, 2.0, 3.0])
    dfs_neg = np.exp(0.02 * times)  # r = -2%

    # Unshifted BK cannot fit negative initial short rate
    model_unshifted = BKTree(0.15, 0.05, 10)
    with pytest.raises(FinError, match="Short rate plus shift must be positive"):
        model_unshifted.build_tree(t_mat, times, dfs_neg)

    # Insufficient shift (shift = 1% for r = -2%) must also raise FinError
    model_small_shift = BKTree(0.15, 0.05, 10, shift=0.01)
    with pytest.raises(FinError, match="Short rate plus shift must be positive"):
        model_small_shift.build_tree(t_mat, times, dfs_neg)

    # Curve starting positive (+1%) then dropping to -3% with shift = 2%
    dfs_later_neg = np.array([1.0, np.exp(-0.01), np.exp(-0.01 + 0.03), np.exp(-0.01 + 0.06)])
    model_mid_fail = BKTree(0.15, 0.05, 6, shift=0.02)
    with pytest.raises(FinError, match="Short rate plus shift must be positive"):
        model_mid_fail.build_tree(t_mat, times, dfs_later_neg)


def test_shifted_bk_parameter_validation_and_repr():
    """Verify input validation for shift and string representation."""
    with pytest.raises(FinError, match="Shift parameter should be float or int"):
        BKTree(0.20, 0.05, 10, shift="invalid")

    model = BKTree(0.20, 0.05, 10, shift=0.03)
    rep = repr(model)
    assert "Black-Karasinski Model" in rep
    assert "shift" in rep
