# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.utils.global_types import OptionTypes
from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.products.equity.equity_cliquet_option import EquityCliquetOption

########################################################################################


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
    discount_curve = DiscountCurveFlat(value_dt, interest_rate)
    dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

    v = cliquet_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(v, 4) == 34.5287
