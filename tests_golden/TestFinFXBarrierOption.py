###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.models.process_simulator import ProcessTypes
from financepy.models.process_simulator import FinGBMNumericalScheme
from financepy.models.black_scholes import BlackScholes
from financepy.products.fx.fx_barrier_option import FinFXBarrierTypes
from financepy.products.fx.fx_barrier_option import FXBarrierOption
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinFXBarrierOption():

    valuation_date = Date(1, 1, 2015)
    expiry_date = Date(1, 1, 2016)
    spot_fx_rate = 100.0
    currency_pair = "USDJPY"
    volatility = 0.20
    dom_interest_rate = 0.05
    forInterestRate = 0.02
    option_type = FinFXBarrierTypes.DOWN_AND_OUT_CALL
    notional = 100.0
    notional_currency = "USD"

    drift = dom_interest_rate - forInterestRate
    scheme = FinGBMNumericalScheme.ANTITHETIC
    process_type = ProcessTypes.GBM
    dom_discount_curve = DiscountCurveFlat(valuation_date, dom_interest_rate)
    for_discount_curve = DiscountCurveFlat(valuation_date, forInterestRate)
    model = BlackScholes(volatility)

    ###########################################################################

    import time
    start = time.time()
    num_observations_per_year = 100

    for option_type in FinFXBarrierTypes:

        testCases.header("Type", "K", "B", "S", "Value",
                         "ValueMC", "TIME", "Diff")

        for spot_fx_rate in range(60, 140, 10):
            B = 110.0
            K = 100.0

            option = FXBarrierOption(expiry_date, K, currency_pair,
                                     option_type, B,
                                     num_observations_per_year,
                                     notional, notional_currency)

            value = option.value(valuation_date, spot_fx_rate,
                                 dom_discount_curve, for_discount_curve, model)

            start = time.time()
            model_params = (spot_fx_rate, drift, volatility, scheme)
            value_mc = option.value_mc(valuation_date, spot_fx_rate,
                                       dom_interest_rate, process_type,
                                       model_params)

            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value

            testCases.print(option_type, K, B, spot_fx_rate, value, value_mc,
                            time_elapsed, diff)

        for spot_fx_rate in range(60, 140, 10):
            B = 100.0
            K = 110.0

            option = FXBarrierOption(expiry_date, K, currency_pair,
                                     option_type, B,
                                     num_observations_per_year,
                                     notional, notional_currency)

            value = option.value(valuation_date, spot_fx_rate,
                                 dom_discount_curve, for_discount_curve, model)

            start = time.time()
            model_params = (spot_fx_rate, drift, volatility, scheme)
            value_mc = option.value_mc(valuation_date,
                                       spot_fx_rate,
                                       dom_interest_rate,
                                       process_type,
                                       model_params)

            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value

            testCases.print(option_type, K, B, spot_fx_rate, value, value_mc,
                            time_elapsed, diff)

    end = time.time()

##########################################################################

    spot_fx_rates = range(50, 150, 50)
    B = 105.0

    testCases.header("Type", "K", "B", "S:", "Value", "Delta", "Vega", "Theta")

    for option_type in FinFXBarrierTypes:
        for spot_fx_rate in spot_fx_rates:
            barrier_option = FXBarrierOption(expiry_date,
                                             100.0,
                                             currency_pair,
                                             option_type,
                                             B,
                                             num_observations_per_year,
                                             notional,
                                             notional_currency)

            value = barrier_option.value(valuation_date,
                                         spot_fx_rate,
                                         dom_discount_curve,
                                         for_discount_curve,
                                         model)

            delta = barrier_option.delta(valuation_date,
                                         spot_fx_rate,
                                         dom_discount_curve,
                                         for_discount_curve,
                                         model)

            vega = barrier_option.vega(valuation_date,
                                       spot_fx_rate,
                                       dom_discount_curve,
                                       for_discount_curve,
                                       model)

            theta = barrier_option.theta(valuation_date,
                                         spot_fx_rate,
                                         dom_discount_curve,
                                         for_discount_curve,
                                         model)

            testCases.print(option_type,
                            K,
                            B,
                            spot_fx_rate,
                            value,
                            delta,
                            vega,
                            theta)

###############################################################################


test_FinFXBarrierOption()
testCases.compareTestCases()
