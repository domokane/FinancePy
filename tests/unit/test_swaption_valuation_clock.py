"""Valuation dates must define a common clock for option and curve times."""

import pytest

from financepy.models.black import Black
from financepy.models.bk_tree import BKTree
from financepy.models.hw_tree import HWTree
from financepy.utils.error import FinError
from financepy.utils.global_vars import G_DAYS_IN_YEAR
from . import test_FinIborSwaption as fixture


def make_option():
    """Use the existing broad swaption fixture with a later settlement date."""
    return fixture.IborSwaption(
        fixture.settle_dt, fixture.exercise_dt, fixture.swap_maturity_dt,
        fixture.SwapTypes.PAY, 0.02, fixture.swap_fixed_freq_type,
        fixture.swap_fixed_day_count_type,
    )


@pytest.mark.parametrize("cash_settled", [False, True])
def test_black_expiry_is_measured_from_valuation_date(monkeypatch, cash_settled):
    """Check both payment conventions with valuation preceding settlement."""
    model = Black(0.2)
    original = model.value
    expiries = []

    def record_value(*args):
        """Record the actual model expiry and preserve its real calculation."""
        expiries.append(args[2])
        return original(*args)

    monkeypatch.setattr(model, "value", record_value)
    option = make_option()
    if cash_settled:
        option.cash_settled_value(fixture.value_dt, fixture.libor_curve, 0.05, model)
    else:
        option.value(fixture.value_dt, fixture.libor_curve, model)

    expected = (fixture.exercise_dt - fixture.value_dt) / G_DAYS_IN_YEAR
    assert expiries == [expected]


def test_tree_horizon_is_measured_from_valuation_date(monkeypatch):
    """Align the tree end time with curve and coupon times."""
    model = BKTree(0.01, 0.01)
    original = model.build_tree
    horizons = []

    def record_tree(*args):
        """Record the horizon and build the actual interest-rate tree."""
        horizons.append(args[0])
        return original(*args)

    monkeypatch.setattr(model, "build_tree", record_tree)
    make_option().value(fixture.value_dt, fixture.libor_curve, model)

    expected = (fixture.swap_maturity_dt - fixture.value_dt) / G_DAYS_IN_YEAR
    assert horizons == [expected]


def test_small_hull_white_volatility_matches_discounted_swap_intrinsic():
    """Check the corrected clock against the deterministic swap cashflows."""
    option = make_option()
    price = option.value(fixture.value_dt, fixture.libor_curve, HWTree(1e-5, 1e-5))
    swap = option.underlying_swap
    expected = max(
        (swap.swap_rate(fixture.value_dt, fixture.libor_curve) - option.fixed_cpn)
        * swap.pv01(fixture.value_dt, fixture.libor_curve) * option.notional,
        0.0,
    ) / fixture.libor_curve.df(fixture.settle_dt)

    assert price == pytest.approx(expected, abs=0.03)
