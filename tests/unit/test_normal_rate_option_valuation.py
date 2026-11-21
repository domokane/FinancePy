"""Independent Gaussian payoff and cash-flow references for normal rate options."""

from math import exp, pi, sqrt

import numpy as np
import pytest
from scipy.integrate import quad

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.bachelier import Bachelier
from financepy.models.black import Black
from financepy.models.sabr_normal import NormalSABR
from financepy.products.rates.ibor_cap_floor import IborCapFloor
from financepy.products.rates.ibor_swaption import IborSwaption
from financepy.utils.calendar import CalendarTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.error import FinError
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.global_types import CapFloorTypes, OptionTypes, SwapTypes


VALUE_DT = Date(1, 1, 2025)
RATE = -0.01
NOTIONAL = 1_000_000.0
ALPHA, RHO, NU = 0.012, -0.40, 0.70


def _years(dt):
    return (dt - VALUE_DT) / 365.0


def _discount(dt):
    return exp(-RATE * _years(dt))


def _normal_sabr_std(forward, strike, expiry, nu):
    # Evaluate x(z)/z as an integral, not with the implementation's asinh
    # difference or small-z Taylor series.
    z = nu * (forward - strike) / ALPHA
    average, _ = quad(
        lambda u: 1.0 / sqrt(1.0 - 2.0 * RHO * z * u + (z * u) ** 2), 0.0, 1.0
    )
    correction = 1.0 + (2.0 - 3.0 * RHO**2) * nu**2 * expiry / 24.0
    return ALPHA * correction / average * sqrt(expiry)


def _normal_payoff(forward, strike, std, sign):
    # Integrate the payoff directly in its exercise region, using symmetry
    # for puts. No calls to the library pricing implementation or CDF.
    intrinsic = sign * (forward - strike)
    if std == 0.0:
        return max(intrinsic, 0.0)
    value, _ = quad(
        lambda x: (intrinsic + std * x) * exp(-x * x / 2.0) / sqrt(2.0 * pi),
        -intrinsic / std,
        np.inf,
        epsabs=1e-14,
    )
    return value


@pytest.mark.parametrize("nu", [0.0, NU])
@pytest.mark.parametrize(
    "forward, strike, expiry",
    [(-0.02, -0.01, 2.0), (-0.01, -0.01, 1.5), (-0.02, 0.01, 0.5),
     (0.04, 0.035, 3.0), (-0.01, -0.012, 0.0)],
)
def test_normal_sabr_value_matches_gaussian_payoff(nu, forward, strike, expiry):
    model = NormalSABR(ALPHA, RHO, nu)
    std = _normal_sabr_std(forward, strike, expiry, nu)
    df = 0.94
    for sign, option in [(1, OptionTypes.EUROPEAN_CALL), (-1, OptionTypes.EUROPEAN_PUT)]:
        expected = df * _normal_payoff(forward, strike, std, sign)
        assert model.value(forward, strike, expiry, df, option) == pytest.approx(
            expected, abs=2e-14
        )


def test_normal_sabr_rejects_negative_approximation():
    model = NormalSABR(ALPHA, 0.99, 3.0)
    with pytest.raises(FinError, match="positive volatility"):
        model.value(-0.01, -0.01, 5.0, 1.0, OptionTypes.EUROPEAN_CALL)


@pytest.mark.filterwarnings("error")
def test_bachelier_broadcasts_at_zero_expiry_without_warning():
    forward = np.array([[-0.02], [-0.01]])
    strike = -0.015
    expiry = np.array([0.0, 1.5])
    val = Bachelier(ALPHA).value(forward, strike, expiry, 0.97, OptionTypes.EUROPEAN_PUT)
    expected = [
        [0.97 * _normal_payoff(f, strike, ALPHA * sqrt(t), -1) for t in expiry]
        for f in forward[:, 0]
    ]
    np.testing.assert_allclose(val, expected, atol=2e-14, rtol=0.0)


@pytest.mark.parametrize("expiry", [VALUE_DT, Date(1, 1, 2027)])
@pytest.mark.parametrize("sabr", [False, True])
def test_physical_swaption_matches_cashflows_and_gaussian_payoff(expiry, sabr):
    model = NormalSABR(ALPHA, RHO, NU) if sabr else Bachelier(ALPHA)
    curve = FlatDiscountCurve(VALUE_DT, RATE)
    settlement = VALUE_DT if expiry == VALUE_DT else Date(1, 2, 2025)
    strike = -0.012
    dates = [Date(1, 1, expiry.y + i) for i in range(4)]
    annuity = sum((d1 - d0) / 365.0 * _discount(d1) for d0, d1 in zip(dates, dates[1:]))
    forward = (_discount(expiry) - _discount(dates[-1])) / annuity
    assert forward < 0.0
    std = _normal_sabr_std(forward, strike, _years(expiry), NU if sabr else 0.0)
    payer_receiver = []
    for sign, side in [(1, SwapTypes.PAY), (-1, SwapTypes.RECEIVE)]:
        swaption = IborSwaption(
            settlement, expiry, dates[-1], side, strike, FrequencyTypes.ANNUAL,
            DayCountTypes.ACT_365F, notional=NOTIONAL, cal_type=CalendarTypes.NONE,
        )
        value = swaption.value(VALUE_DT, curve, model)
        multiplier = NOTIONAL * annuity / _discount(settlement)
        # FlatDiscountCurve floors the anchor time at G_SMALL; allow its
        # roughly 1e-8 cash impact on this notional.
        assert value == pytest.approx(
            multiplier * _normal_payoff(forward, strike, std, sign), abs=1e-7
        )
        payer_receiver.append(value)
    assert payer_receiver[0] - payer_receiver[1] == pytest.approx(
        NOTIONAL * (_discount(expiry) - _discount(dates[-1]) - strike * annuity)
        / _discount(settlement), abs=1e-7
    )


@pytest.mark.parametrize("strike", [-0.012, 0.0])
@pytest.mark.parametrize("sabr", [False, True])
def test_forward_cap_floor_matches_sum_of_caplet_payoffs(strike, sabr):
    model = NormalSABR(ALPHA, RHO, NU) if sabr else Bachelier(ALPHA)
    curve = FlatDiscountCurve(VALUE_DT, RATE)
    dates = [Date(1, 1, 2026), Date(1, 4, 2026), Date(1, 7, 2026)]
    for sign, side in [(1, CapFloorTypes.CAP), (-1, CapFloorTypes.FLOOR)]:
        cap_floor = IborCapFloor(
            dates[0], dates[-1], side, strike, accrual_dc_type=DayCountTypes.ACT_365F,
            notional=NOTIONAL, cal_type=CalendarTypes.NONE,
        )
        expected_legs = []
        for start, end in zip(dates, dates[1:]):
            accrual = (end - start) / 365.0
            forward = (_discount(start) / _discount(end) - 1.0) / accrual
            std = _normal_sabr_std(forward, strike, _years(start), NU if sabr else 0.0)
            expected_legs.append(
                _discount(end) * accrual * NOTIONAL * _normal_payoff(forward, strike, std, sign)
            )
        total = cap_floor.value(VALUE_DT, curve, model)
        assert total == pytest.approx(sum(expected_legs), abs=1e-8)
        np.testing.assert_allclose(
            cap_floor.caplet_floorlet_values[1:], expected_legs, atol=1e-8, rtol=0.0
        )


def test_black_cap_floor_still_rejects_negative_strike():
    cap = IborCapFloor(Date(1, 1, 2026), Date(1, 7, 2026), CapFloorTypes.CAP, -0.01)
    with pytest.raises(FinError, match="Strike < 0.0"):
        cap.value(VALUE_DT, FlatDiscountCurve(VALUE_DT, RATE), Black(0.2))
