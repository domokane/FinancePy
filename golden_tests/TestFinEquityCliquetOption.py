###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.products.equity.equity_cliquet_option import EquityCliquetOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes
from FinTestCases import FinTestCases, globalTestCaseMode

test_cases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_EquityCliquetOption():

    start_dt = Date(1, 1, 2014)
    final_expiry_dt = Date(1, 1, 2017)
    freq_type = FrequencyTypes.QUARTERLY
    option_type = OptionTypes.EUROPEAN_CALL

    cliquetOption = EquityCliquetOption(start_dt,
                                        final_expiry_dt,
                                        option_type,
                                        freq_type)

    value_dt = Date(1, 1, 2015)
    stock_price = 100.0
    volatility = 0.20
    interest_rate = 0.05
    dividend_yield = 0.02
    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(value_dt, interest_rate)
    dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

    v = cliquetOption.value(value_dt,
                            stock_price,
                            discount_curve,
                            dividend_curve,
                            model)

    test_cases.header("LABEL", "VALUE")
    test_cases.print("FINANCEPY", v)

###############################################################################


test_EquityCliquetOption()
test_cases.compareTestCases()
