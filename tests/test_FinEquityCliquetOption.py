###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.global_types import OptionTypes
from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.products.equity.equity_cliquet_option import EquityCliquetOption


def test_EquityCliquetOption():

    start_date = Date(1, 1, 2014)
    final_expiry_date = Date(1, 1, 2017)
    freq_type = FrequencyTypes.QUARTERLY
    option_type = OptionTypes.EUROPEAN_CALL

    cliquetOption = EquityCliquetOption(start_date,
                                        final_expiry_date,
                                        option_type,
                                        freq_type)

    valuation_date = Date(1, 1, 2015)
    stock_price = 100.0
    volatility = 0.20
    interest_rate = 0.05
    dividend_yield = 0.02
    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

    v = cliquetOption.value(valuation_date,
                            stock_price,
                            discount_curve,
                            dividend_curve,
                            model)

    assert round(v, 4) == 34.5287
