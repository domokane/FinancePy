import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.products.equity.equity_fixed_lookback_option import EquityFixedLookbackOption
from financepy.products.equity.equity_float_lookback_option import EquityFloatLookbackOption
from financepy.products.fx.fx_fixed_lookback_option import FXFixedLookbackOption
from financepy.products.fx.fx_float_lookback_option import FXFloatLookbackOption
from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes

# With w = 2(r - q) / vol^2 beyond +-100 the term (S/X)^(-w) N(.) is still of order one when
# ln(X/S) is close to (r - q) T, so it cannot be dropped. Reference values are from QuantLib's
# AnalyticContinuousFixedLookbackEngine and AnalyticContinuousFloatingLookbackEngine.

VALUE_DT = Date(1, 1, 2025)
EXPIRY_DT = Date(1, 7, 2025)
STOCK_PRICE = 100.0

# (option type, strike, r, q, vol, running extreme, reference value)
FIXED_CASES = [
    (OptionTypes.EUROPEAN_CALL, 105.0, 0.10, 0.0, 0.04, 100.0, 1.491708),  # w = 125, strike above max
    (OptionTypes.EUROPEAN_CALL, 95.0, 0.10, 0.0, 0.04, 103.0, 10.520411),  # w = 125, strike below max
    (OptionTypes.EUROPEAN_PUT, 95.0, 0.0, 0.08, 0.03, 100.0, 0.483659),  # w = -178, strike below min
    (OptionTypes.EUROPEAN_PUT, 100.0, 0.0, 0.08, 0.03, 97.0, 4.639173),  # w = -178, strike above min
]

# floating put, r = 0.10, q = 0, vol = 0.04 (w = 125), running maximum 105
FLOAT_PUT_REFERENCE = 1.411852


@pytest.mark.parametrize("opt_type, strike, r, q, vol, extreme, expected", FIXED_CASES)
def test_equity_fixed_lookback_large_carry_to_variance(opt_type, strike, r, q, vol, extreme, expected):
    option = EquityFixedLookbackOption(EXPIRY_DT, opt_type, strike)
    discount = FlatDiscountCurve(VALUE_DT, r)
    dividend = FlatDiscountCurve(VALUE_DT, q)

    price = option.value(VALUE_DT, STOCK_PRICE, discount, dividend, BlackScholes(vol), extreme)

    assert price == pytest.approx(expected, abs=1.0e-5)


@pytest.mark.parametrize("opt_type, strike, r, q, vol, extreme, expected", FIXED_CASES)
def test_fx_fixed_lookback_large_carry_to_variance(opt_type, strike, r, q, vol, extreme, expected):
    option = FXFixedLookbackOption(EXPIRY_DT, opt_type, strike)
    domestic = FlatDiscountCurve(VALUE_DT, r)
    foreign = FlatDiscountCurve(VALUE_DT, q)

    price = option.value(VALUE_DT, STOCK_PRICE, domestic, foreign, vol, extreme)

    assert price == pytest.approx(expected, abs=1.0e-5)


def test_equity_floating_lookback_put_large_carry_to_variance():
    option = EquityFloatLookbackOption(EXPIRY_DT, OptionTypes.EUROPEAN_PUT)
    discount = FlatDiscountCurve(VALUE_DT, 0.10)
    dividend = FlatDiscountCurve(VALUE_DT, 0.0)

    price = option.value(VALUE_DT, STOCK_PRICE, discount, dividend, BlackScholes(0.04), 105.0)

    assert price == pytest.approx(FLOAT_PUT_REFERENCE, abs=1.0e-5)


def test_fx_floating_lookback_put_large_carry_to_variance():
    option = FXFloatLookbackOption(EXPIRY_DT, OptionTypes.EUROPEAN_PUT)
    domestic = FlatDiscountCurve(VALUE_DT, 0.10)
    foreign = FlatDiscountCurve(VALUE_DT, 0.0)

    price = option.value(VALUE_DT, STOCK_PRICE, domestic, foreign, 0.04, 105.0)

    assert price == pytest.approx(FLOAT_PUT_REFERENCE, abs=1.0e-5)
