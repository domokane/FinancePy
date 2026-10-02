# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
import pytest
from financepy.utils.error import FinError
from financepy.utils.global_types import OptionTypes
from financepy.models.black_scholes_analytic import european_value

from financepy.utils.date import Date
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.products.equity.equity_chooser_option import EquityChooserOption


@pytest.mark.parametrize("value_dt", ["2026-01-01", Date(2, 4, 2026)])
def test_mc_rejects_invalid_or_late_valuation_dates(value_dt):
    """Invalid dates cannot reach negative Monte Carlo time intervals."""
    curve = FlatDiscountCurve(Date(1, 1, 2026), 0.03)
    option = EquityChooserOption(Date(1, 4, 2026), Date(1, 7, 2026), Date(1, 8, 2026), 100.0, 100.0)
    with pytest.raises(FinError):
        option.value_mc(value_dt, 100.0, curve, curve, BlackScholes(0.2))


@pytest.mark.parametrize("wrong_curve", ["discount", "dividend"])
def test_mc_rejects_curve_date_mismatch(wrong_curve):
    """Both simulation curves must describe the requested valuation date."""
    value_dt = Date(1, 1, 2026)
    curve = FlatDiscountCurve(value_dt, 0.03)
    wrong = FlatDiscountCurve(value_dt.add_days(1), 0.03)
    option = EquityChooserOption(Date(1, 4, 2026), Date(1, 7, 2026), Date(1, 8, 2026), 100.0, 100.0)
    discount = wrong if wrong_curve == "discount" else curve
    dividend = wrong if wrong_curve == "dividend" else curve
    with pytest.raises(FinError):
        option.value_mc(value_dt, 100.0, discount, dividend, BlackScholes(0.2))


def test_mc_on_choose_date_matches_independent_best_vanilla_option():
    """At the decision date the holder chooses the more valuable vanilla."""
    choose_dt = Date(1, 4, 2026)
    call_expiry = Date(1, 7, 2026)
    put_expiry = Date(1, 8, 2026)
    option = EquityChooserOption(choose_dt, call_expiry, put_expiry, 100.0, 100.0)
    discount = FlatDiscountCurve(choose_dt, 0.03)
    dividend = FlatDiscountCurve(choose_dt, 0.01)
    call = european_value(100.0, (call_expiry - choose_dt) / 365.0, 100.0, .03, .01, .2, OptionTypes.EUROPEAN_CALL.value)
    put = european_value(100.0, (put_expiry - choose_dt) / 365.0, 100.0, .03, .01, .2, OptionTypes.EUROPEAN_PUT.value)
    assert option.value_mc(choose_dt, 100.0, discount, dividend, BlackScholes(.2)) == pytest.approx(max(call, put))


def assert_close(value, expected, tol=2.0e-3):
    assert np.isclose(value, expected, atol=tol), (
        f"value={value:.10f}, expected={expected:.10f}, " f"diff={value - expected:.10f}"
    )


########################################################################################


def test_equity_chooser_option_haug():
    """Following example in Haug Page 130"""

    value_dt = Date(1, 1, 2015)
    choose_dt = Date(2, 4, 2015)
    call_expiry_dt = Date(1, 7, 2015)
    put_expiry_dt = Date(2, 8, 2015)
    call_strike = 55.0
    put_strike = 48.0
    stock_price = 50.0
    volatility = 0.35
    interest_rate = 0.10
    dividend_yield = 0.05

    model = BlackScholes(volatility)
    discount_curve = FlatDiscountCurve(value_dt, interest_rate)
    dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

    chooser_option = EquityChooserOption(choose_dt, call_expiry_dt, put_expiry_dt, call_strike, put_strike)

    v = chooser_option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    v_mc = chooser_option.value_mc(value_dt, stock_price, discount_curve, dividend_curve, model, 20000)

    v_haug = 6.0508

    assert_close(v, 6.034)
    assert_close(v_haug, 6.0508)
    assert_close(v_mc, 6.0320)


########################################################################################


def test_equity_chooser_option_matlab():
    """https://fr.mathworks.com/help/fininst/chooserbybls.html"""

    value_dt = Date(1, 6, 2007)
    choose_date = Date(31, 8, 2007)
    call_expiry_dt = Date(2, 12, 2007)
    put_expiry_dt = Date(2, 12, 2007)
    call_strike = 60.0
    put_strike = 60.0
    stock_price = 50.0
    volatility = 0.20
    interest_rate = 0.10
    dividend_yield = 0.05

    model = BlackScholes(volatility)

    discount_curve = FlatDiscountCurve(value_dt, interest_rate)
    dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

    chooser_option = EquityChooserOption(choose_date, call_expiry_dt, put_expiry_dt, call_strike, put_strike)

    v = chooser_option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    v_mc = chooser_option.value_mc(value_dt, stock_price, discount_curve, dividend_curve, model, 20000)

    v_matlab = 8.9308

    assert_close(v, 8.931)
    assert_close(v_matlab, 8.931)
    assert_close(v_mc, 8.927)


########################################################################################


def test_equity_chooser_option_derivicom():
    """http://derivicom.com/support/finoptionsxl/index.html?complex_chooser.htm"""

    value_dt = Date(1, 1, 2007)
    choose_date = Date(1, 2, 2007)
    call_expiry_dt = Date(1, 4, 2007)
    put_expiry_dt = Date(1, 5, 2007)
    call_strike = 40.0
    put_strike = 35.0
    stock_price = 38.0
    volatility = 0.20
    interest_rate = 0.08
    dividend_yield = 0.0625

    model = BlackScholes(volatility)
    discount_curve = FlatDiscountCurve(value_dt, interest_rate)
    dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

    chooser_option = EquityChooserOption(choose_date, call_expiry_dt, put_expiry_dt, call_strike, put_strike)

    v = chooser_option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    v_mc = chooser_option.value_mc(value_dt, stock_price, discount_curve, dividend_curve, model, 20000)

    v_derivicom = 1.0989

    assert_close(v, 1.105)
    assert_close(v_derivicom, 1.0989)
    assert_close(v_mc, 1.1046)
