########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

import numpy as np
import pytest

from financepy.utils.global_types import OptionTypes
from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.models.black_scholes_analytic import european_value
from financepy.products.equity.equity_forward_start_option import (
    EquityForwardStartOption,
)
from financepy.products.equity.equity_cliquet_option import EquityCliquetOption

value_dt = Date(1, 1, 2026)
start_dt = Date(1, 7, 2026)
expiry_dt = Date(1, 7, 2027)
stock_price = 100.0
volatility = 0.25
interest_rate = 0.05
dividend_yield = 0.02
model = BlackScholes(volatility)
discount_curve = FlatDiscountCurve(value_dt, interest_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)


########################################################################################


def _rubinstein(opt_type, t_start, t_expiry):
    """Rubinstein (1991) at-the-money forward-start option with continuous rates."""
    tau = t_expiry - t_start
    v = european_value(1.0, tau, 1.0, interest_rate, dividend_yield, volatility, opt_type.value)
    return stock_price * np.exp(-dividend_yield * t_start) * v


def test_forward_start_matches_rubinstein():
    """A single period reproduces the closed form for a call and a put."""
    t_start = (start_dt - value_dt) / 365.0
    t_expiry = (expiry_dt - value_dt) / 365.0
    for opt_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
        option = EquityForwardStartOption(start_dt, expiry_dt, opt_type, FrequencyTypes.ANNUAL)
        v = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)
        assert v == pytest.approx(_rubinstein(opt_type, t_start, t_expiry), rel=1e-10)
        assert len(option.v_options) == 1
        assert option.actual_dts == [expiry_dt]


def test_forward_start_is_a_cliquet_from_the_start_date():
    """With several periods the value is the sum of the forward-start options, and
    equals the cliquet option valued on the start date deflated by the dividend
    yield to that date."""
    option = EquityForwardStartOption(start_dt, expiry_dt, OptionTypes.EUROPEAN_CALL, FrequencyTypes.QUARTERLY)
    v = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)
    assert len(option.v_options) == 4
    assert v == pytest.approx(sum(option.v_options))

    cliquet = EquityCliquetOption(start_dt, expiry_dt, OptionTypes.EUROPEAN_CALL, FrequencyTypes.QUARTERLY)
    v_cliquet = cliquet.value(
        start_dt,
        stock_price,
        FlatDiscountCurve(start_dt, interest_rate),
        FlatDiscountCurve(start_dt, dividend_yield),
        model,
    )
    t_start = (start_dt - value_dt) / 365.0
    assert v == pytest.approx(v_cliquet * np.exp(-dividend_yield * t_start), rel=1e-10)


def test_forward_start_rejects_valuation_after_start():
    option = EquityForwardStartOption(start_dt, expiry_dt, OptionTypes.EUROPEAN_CALL, FrequencyTypes.ANNUAL)
    with pytest.raises(Exception, match="strike has been set"):
        option.value(start_dt.add_days(1), stock_price, FlatDiscountCurve(start_dt.add_days(1), interest_rate), FlatDiscountCurve(start_dt.add_days(1), dividend_yield), model)
    with pytest.raises(Exception, match="after start date"):
        EquityForwardStartOption(expiry_dt, start_dt, OptionTypes.EUROPEAN_CALL, FrequencyTypes.ANNUAL)


def test_forward_start_repr():
    option = EquityForwardStartOption(start_dt, expiry_dt, OptionTypes.EUROPEAN_PUT, FrequencyTypes.ANNUAL)
    text = repr(option)
    assert "EquityForwardStartOption" in text
    assert "CALENDAR TYPE" in text
