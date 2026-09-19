# Copyright (C) 2018, 2019, 2020 Dominic O'Kane



# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path
_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause
_install_double_click_pause()
import add_fp_to_path

from financepy.products.equity.equity_compound_option import EquityCompoundOption
from financepy.utils.global_types import OptionTypes
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.date import Date


########################################################################################


def test_equity_compound_option():

    value_dt = Date(1, 1, 2015)
    expiry_dt1 = Date(1, 1, 2017)
    expiry_dt2 = Date(1, 1, 2018)
    k1 = 5.0
    k2 = 95.0
    stock_price = 85.0
    volatility = 0.15
    interest_rate = 0.035
    dividend_yield = 0.01

    model = BlackScholes(volatility)
    discount_curve = FlatDiscountCurve(value_dt, interest_rate)
    dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

    num_steps_list = [100, 200, 500, 1000]

    stock_price = 85.0

    print(
        "TYPE1", "TYPE2", "k_1", "k_2", "S", "TreeSteps", "Exact", "TreeValue"
    )

    for opt_type1 in [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT]:
        for opt_type2 in [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT]:

            cmpd_option = EquityCompoundOption(
                expiry_dt1, opt_type1, k1, expiry_dt2, opt_type2, k2
            )

            for num_steps in num_steps_list:

                value = cmpd_option.value(
                    value_dt,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model,
                )

                values_tree = cmpd_option.value_tree(value_dt,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model,
                    num_steps,
                )

                print(
                    opt_type1,
                    opt_type2,
                    k1,
                    k2,
                    stock_price,
                    num_steps,
                    value,
                    values_tree[0],
                )

    stock_price = 85.0

    print(
        "TYPE1", "TYPE2", "k_1", "k_2", "S", "TreeSteps", "Exact", "TreeValue"
    )

    for opt_type1 in [OptionTypes.AMERICAN_CALL, OptionTypes.AMERICAN_PUT]:
        for opt_type2 in [OptionTypes.AMERICAN_CALL, OptionTypes.AMERICAN_PUT]:

            cmpd_option = EquityCompoundOption(
                expiry_dt1, opt_type1, k1, expiry_dt2, opt_type2, k2
            )

            for num_steps in num_steps_list:

                value = cmpd_option.value(
                    value_dt,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model,
                    num_steps,
                )

                values_tree = cmpd_option.value_tree(
                    value_dt,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model,
                    num_steps,
                )

                print(
                    opt_type1,
                    opt_type2,
                    k1,
                    k2,
                    stock_price,
                    num_steps,
                    value,
                    values_tree[0],
                )

    print(
        "TYPE1",
        "TYPE2",
        "k_1",
        "k_2",
        "S",
        "Exact",
        "TreeSteps",
        "TreeValue",
        "Diff",
        "DELTA",
        "GAMMA",
        "THETA",
    )

    for opt_type1 in [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT]:
        for opt_type2 in [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT]:

            cmpd_option = EquityCompoundOption(
                expiry_dt1, opt_type1, k1, expiry_dt2, opt_type2, k2
            )
            stock_prices = range(70, 100, 10)

            for stock_price in stock_prices:
                value = cmpd_option.value(
                    value_dt,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model,
                )
                delta = cmpd_option.delta(
                    value_dt,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model,
                )
                vega = cmpd_option.vega(
                    value_dt,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model,
                )
                theta = cmpd_option.theta(
                    value_dt,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model,
                )

                values = cmpd_option.value_tree(
                    value_dt,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model,
                )

                diff = value - values[0]

                print(
                    opt_type1,
                    opt_type2,
                    k1,
                    k2,
                    stock_price,
                    value,
                    num_steps,
                    values[0],
                    diff,
                    delta,
                    vega,
                    theta,
                )


########################################################################################

test_equity_compound_option()
