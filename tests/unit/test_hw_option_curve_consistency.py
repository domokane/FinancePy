"""Hull-White rate option PVs must use the supplied curve at cash-flow dates."""

from math import exp, expm1, pi, sqrt

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq

from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.interpolator import InterpTypes
from financepy.models.hw_tree import HWTree
from financepy.products.rates.ibor_cap_floor import IborCapFloor
from financepy.products.rates.ibor_swaption import IborSwaption
from financepy.utils.calendar import CalendarTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.global_types import CapFloorTypes, SwapTypes


VALUE_DT = Date(1, 1, 2025)
EXERCISE_DT = Date(1, 7, 2027)
NOTIONAL = 1_000_000.0
A = 0.07
CURVE_TYPES = [
    InterpTypes.LINEAR_ZERO_RATES,
    InterpTypes.MONOTONE_CONVEX,
    InterpTypes.FLAT_FWD_RATES,
]


def _time(date):
    return (date - VALUE_DT) / 365.0


def _curve(interp_type):
    knots = [VALUE_DT, Date(1, 1, 2026), Date(1, 1, 2029), Date(1, 1, 2033)]
    times = np.array([_time(date) for date in knots])
    dfs = np.exp(-0.01 * times - 0.004 * times**2)
    return DiscountCurve(VALUE_DT, knots, dfs, interp_type)


def _expiry_measure_payoffs(curve, expiry_dt, flow_dates, amounts, strike, sigma, a=A):
    # Under the expiry forward measure, every P(T,U) is a lognormal
    # multiple of the same N(0,1), and E[P(T,U)]=P(0,U)/P(0,T).
    expiry_t = _time(expiry_dt)
    expiry_df = curve.df(expiry_dt)
    std_x = sigma * sqrt(-expm1(-2 * a * expiry_t) / (2 * a) if a else expiry_t)
    dfs = np.array([curve.df(d) / expiry_df for d in flow_dates])
    dt = np.array([_time(d) - expiry_t for d in flow_dates])
    std = (-np.expm1(-a * dt) / a if a else dt) * std_x
    weights = np.asarray(amounts) * dfs * np.exp(-0.5 * std**2)

    def payoff(z):
        return np.dot(weights, np.exp(-std * z)) - strike

    if std_x == 0.0:
        deterministic = payoff(0.0)
        return expiry_df * max(deterministic, 0.0), expiry_df * max(-deterministic, 0.0)

    # Integrate max(payoff, 0) directly, not a sum of ZCB Black prices.
    lo, hi = -12.0, 12.0
    root = brentq(payoff, lo, hi)

    def integrand(z):
        return payoff(z) * exp(-z * z / 2.0) / sqrt(2.0 * pi)

    positive = quad(integrand, lo, root, epsabs=1e-13)[0]
    negative = -quad(integrand, root, hi, epsabs=1e-13)[0]
    return expiry_df * positive, expiry_df * negative


@pytest.mark.parametrize("interp_type", CURVE_TYPES)
@pytest.mark.parametrize("a", [A, 1e-11, 0.0])
@pytest.mark.parametrize("sigma", [0.0, 0.018])
def test_hw_swaption_respects_off_knot_discounts(interp_type, sigma, a):
    curve = _curve(interp_type)
    model = HWTree(sigma, a)
    dates = [Date(1, 7, 2027 + i) for i in range(4)]
    rate = 0.04
    amounts = np.array([rate * (d1 - d0) / 365.0 for d0, d1 in zip(dates, dates[1:])])
    amounts[-1] += 1.0
    call, put = _expiry_measure_payoffs(curve, dates[0], dates[1:], amounts, 1.0, sigma, a)
    for side, expected in [(SwapTypes.RECEIVE, call), (SwapTypes.PAY, put)]:
        product = IborSwaption(
            VALUE_DT, dates[0], dates[-1], side, rate, FrequencyTypes.ANNUAL,
            DayCountTypes.ACT_365F, notional=NOTIONAL, cal_type=CalendarTypes.NONE,
        )
        assert product.value(VALUE_DT, curve, model) == pytest.approx(
            NOTIONAL * expected, abs=1.0
        )


@pytest.mark.parametrize("interp_type", CURVE_TYPES)
@pytest.mark.parametrize("a", [A, 1e-11, 0.0])
@pytest.mark.parametrize("sigma", [0.0, 0.018])
def test_hw_caplet_respects_off_knot_discounts(interp_type, sigma, a):
    curve = _curve(interp_type)
    model = HWTree(sigma, a)
    start, end = EXERCISE_DT, Date(1, 10, 2027)
    rate = 0.033
    scale = 1.0 + rate * (end - start) / 365.0
    call, put = _expiry_measure_payoffs(curve, start, [end], [1.0], 1.0 / scale, sigma, a)
    for side, expected in [(CapFloorTypes.CAP, put), (CapFloorTypes.FLOOR, call)]:
        cap_floor = IborCapFloor(
            start, end, side, rate, accrual_dc_type=DayCountTypes.ACT_365F,
            notional=NOTIONAL, cal_type=CalendarTypes.NONE,
        )
        assert cap_floor.value(VALUE_DT, curve, model) == pytest.approx(
            NOTIONAL * scale * expected, abs=0.15
        )
