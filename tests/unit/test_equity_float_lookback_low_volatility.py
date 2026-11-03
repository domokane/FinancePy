import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.products.equity.equity_float_lookback_option import EquityFloatLookbackOption
from financepy.products.fx.fx_float_lookback_option import FXFloatLookbackOption
from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes


def test_equity_floating_lookback_call_avoids_intermediate_overflow():
    value_dt = Date(1, 1, 2026)
    expiry_dt = Date(1, 1, 2027)
    discount = FlatDiscountCurve(value_dt, 0.02)
    dividend = FlatDiscountCurve(value_dt, 0.07)
    option = EquityFloatLookbackOption(expiry_dt, OptionTypes.EUROPEAN_CALL)

    price = option.value(value_dt, 100.0, discount, dividend, BlackScholes(0.01), 40.0)

    assert price == pytest.approx(54.031435058324617, abs=2.0e-6)


def test_fx_floating_lookback_call_avoids_intermediate_overflow():
    value_dt = Date(1, 1, 2026)
    expiry_dt = Date(1, 1, 2027)
    domestic = FlatDiscountCurve(value_dt, 0.02)
    foreign = FlatDiscountCurve(value_dt, 0.07)
    option = FXFloatLookbackOption(expiry_dt, OptionTypes.EUROPEAN_CALL)

    price = option.value(value_dt, 100.0, domestic, foreign, 0.01, 40.0)

    assert price == pytest.approx(54.031435058324617, abs=2.0e-6)
