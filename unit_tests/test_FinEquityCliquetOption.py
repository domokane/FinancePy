# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

from financepy.utils.global_types import OptionTypes
from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.products.equity.equity_cliquet_option import EquityCliquetOption

########################################################################################


def assert_close(value, expected, tol=2.0e-3):
    assert np.isclose(value, expected, atol=tol), (
        f"value={value:.10f}, expected={expected:.10f}, " f"diff={value - expected:.10f}"
    )


def test_equity_cliquet_option():

    start_dt = Date(1, 1, 2014)
    final_expiry_dt = Date(1, 1, 2017)
    freq_type = FrequencyTypes.QUARTERLY
    opt_type = OptionTypes.EUROPEAN_CALL

    cliquet_option = EquityCliquetOption(start_dt, final_expiry_dt, opt_type, freq_type)

    value_dt = Date(1, 1, 2015)
    stock_price = 100.0
    volatility = 0.20
    interest_rate = 0.05
    dividend_yield = 0.02
    model = BlackScholes(volatility)
    discount_curve = FlatDiscountCurve(value_dt, interest_rate)
    dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

    v = cliquet_option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    assert_close(v, 34.531)
