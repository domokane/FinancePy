###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
from financepy.utils.global_types import OptionTypes
from financepy.products.fx.fx_vanilla_option import FXVanillaOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.day_count import DayCountTypes
from financepy.utils.calendar import CalendarTypes
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.utils.date import Date
import sys
sys.path.append("./..")


def test_FinFXVanillaOptionWystupExample1():
    # Example from Book extract by Uwe Wystup with results in Table 1.2
    # https://mathfinance.com/wp-content/uploads/2017/06/FXOptionsStructuredProducts2e-Extract.pdf

    # Not exactly T=1.0 but close so don't exact exact agreement
    # (in fact I do not get exact agreement even if I do set T=1.0)
    valuation_date = Date(13, 2, 2018)
    expiry_date = Date(13, 2, 2019)

    # In BS the FX rate is the price in domestic of one unit of foreign
    # In case of EURUSD = 1.3 the domestic currency is USD and foreign is EUR
    # DOM = USD , FOR = EUR
    ccy1 = "EUR"
    ccy2 = "USD"
    ccy1CCRate = 0.030  # EUR
    ccy2CCRate = 0.025  # USD

    currency_pair = ccy1 + ccy2  # Always ccy1ccy2
    spot_fx_rate = 1.20
    strike_fx_rate = 1.250
    volatility = 0.10

    notional = 1000000.0

    dom_discount_curve = DiscountCurveFlat(valuation_date, ccy2CCRate)
    for_discount_curve = DiscountCurveFlat(valuation_date, ccy1CCRate)

    model = BlackScholes(volatility)

    # Two examples to show that changing the notional currency and notional
    # keeps the value unchanged
    notional = 1000000.0
    call_option = FXVanillaOption(expiry_date,
                                  strike_fx_rate,
                                  currency_pair,
                                  OptionTypes.EUROPEAN_CALL,
                                  notional,
                                  "EUR", 2)

    value = call_option.value(valuation_date,
                              spot_fx_rate,
                              dom_discount_curve,
                              for_discount_curve,
                              model)

    notional = 1250000.0
    call_option = FXVanillaOption(expiry_date,
                                  strike_fx_rate,
                                  currency_pair,
                                  OptionTypes.EUROPEAN_CALL,
                                  notional,
                                  "USD", 2)

    value = call_option.value(
        valuation_date,
        spot_fx_rate,
        dom_discount_curve,
        for_discount_curve,
        model)

    assert round(value['v'], 4) == 0.0251
    assert round(value['cash_dom'], 4) == 25125.1772
    assert round(value['cash_for'], 4) == 20937.6477
    assert round(value['pips_dom'], 4) == 0.0251
    assert round(value['pips_for'], 4) == 0.0168
    assert round(value['pct_dom'], 4) == 0.0201
    assert round(value['pct_for'], 4) == 0.0209
    assert round(value['not_dom'], 4) == 1250000.0
    assert round(value['not_for'], 4) == 1000000.0
    assert value['ccy_dom'] == 'USD'
    assert value['ccy_for'] == 'EUR'

    delta = call_option.delta(
        valuation_date,
        spot_fx_rate,
        dom_discount_curve,
        for_discount_curve,
        model)

    assert round(delta['pips_spot_delta'], 4) == 0.3315
    assert round(delta['pips_fwd_delta'], 4) == 0.3416
    assert round(delta['pct_spot_delta_prem_adj'], 4) == 0.3105
    assert round(delta['pct_fwd_delta_prem_adj'], 4) == 0.3200


###############################################################################


def test_FinFXVanillaOptionWystupExample2():
    # Example Bloomberg Pricing at
    # https://stackoverflow.com/questions/48778712/fx-vanilla-call-price-in-quantlib-doesnt-match-bloomberg

    valuation_date = Date(13, 2, 2018)
    expiry_date = Date(13, 2, 2019)

    # In BS the FX rate is the price in domestic of one unit of foreign
    # In case of EURUSD = 1.3 the domestic currency is USD and foreign is EUR
    # DOM = USD , FOR = EUR
    ccy1 = "EUR"
    ccy2 = "USD"
    ccy1CCRate = 0.0396  # EUR
    ccy2CCRate = 0.0357  # USD

    currency_pair = ccy1 + ccy2  # Always ccy1ccy2
    spot_fx_rate = 0.9090
    strike_fx_rate = 0.9090
    volatility = 0.12

    notional = 1000000.0

    dom_discount_curve = DiscountCurveFlat(valuation_date, ccy2CCRate)
    for_discount_curve = DiscountCurveFlat(valuation_date, ccy1CCRate)

    model = BlackScholes(volatility)

    # Two examples to show that changing the notional currency and notional
    # keeps the value unchanged
    notional = 1000000.0
    call_option = FXVanillaOption(expiry_date,
                                  strike_fx_rate,
                                  currency_pair,
                                  OptionTypes.EUROPEAN_PUT,
                                  notional,
                                  "EUR", 2)

    value = call_option.value(
        valuation_date,
        spot_fx_rate,
        dom_discount_curve,
        for_discount_curve,
        model)

    assert round(value['v'], 4) == 0.0436
    assert round(value['cash_dom'], 4) == 43612.8769
    assert round(value['cash_for'], 4) == 47978.9625
    assert round(value['pips_dom'], 4) == 0.0436
    assert round(value['pips_for'], 4) == 0.0528
    assert round(value['pct_dom'], 4) == 0.0480
    assert round(value['pct_for'], 4) == 0.0480
    assert round(value['not_dom'], 4) == 909000.0
    assert round(value['not_for'], 4) == 1000000.0
    assert value['ccy_dom'] == 'USD'
    assert value['ccy_for'] == 'EUR'

    delta = call_option.delta(
        valuation_date,
        spot_fx_rate,
        dom_discount_curve,
        for_discount_curve,
        model)

    assert round(delta['pips_spot_delta'], 4) == -0.4700
    assert round(delta['pips_fwd_delta'], 4) == -0.4890
    assert round(delta['pct_spot_delta_prem_adj'], 4) == -0.5180
    assert round(delta['pct_fwd_delta_prem_adj'], 4) == -0.5389

###############################################################################


def test_FinFXVanillaOptionBloombergExample():
    # Example Bloomberg Pricing at
    # https://stackoverflow.com/questions/48778712/fx-vanilla-call-price-in-quantlib-doesnt-match-bloomberg

    valuation_date = Date(13, 2, 2018)
    expiry_date = Date(15, 2, 2019)

    # In BS the FX rate is the price in domestic of one unit of foreign
    # In case of EURUSD = 1.3 the domestic currency is USD and foreign is EUR
    # DOM = USD , FOR = EUR
    forName = "EUR"
    domName = "USD"
    forDepoRate = 0.05  # EUR
    domDepoRate = 0.02  # USD

    currency_pair = forName + domName  # Always FORDOM
    spot_fx_rate = 1.30
    strike_fx_rate = 1.3650
    volatility = 0.20

    spot_days = 0
    settlement_date = valuation_date.add_weekdays(spot_days)
    maturity_date = settlement_date.add_months(12)
    notional = 1000000.0
    notional_currency = "EUR"
    calendar_type = CalendarTypes.TARGET

    depos = []
    fras = []
    swaps = []
    depo = IborDeposit(settlement_date, maturity_date, domDepoRate,
                       DayCountTypes.ACT_360, notional, calendar_type)
    depos.append(depo)
    dom_discount_curve = IborSingleCurve(valuation_date, depos, fras, swaps)

    depos = []
    fras = []
    swaps = []
    depo = IborDeposit(settlement_date, maturity_date, forDepoRate,
                       DayCountTypes.ACT_360, notional, calendar_type)
    depos.append(depo)
    for_discount_curve = IborSingleCurve(valuation_date, depos, fras, swaps)

    model = BlackScholes(volatility)

    call_option = FXVanillaOption(expiry_date,
                                  strike_fx_rate,
                                  currency_pair,
                                  OptionTypes.EUROPEAN_CALL,
                                  notional,
                                  notional_currency, 2)

    value = call_option.value(
        valuation_date,
        spot_fx_rate,
        dom_discount_curve,
        for_discount_curve,
        model)

    assert round(value['v'], 4) == 0.0601
    assert round(value['cash_dom'], 4) == 60145.5078
    assert round(value['cash_for'], 4) == 46265.7752
    assert round(value['pips_dom'], 4) == 0.0601
    assert round(value['pips_for'], 4) == 0.0339
    assert round(value['pct_dom'], 4) == 0.0441
    assert round(value['pct_for'], 4) == 0.0463
    assert round(value['not_dom'], 4) == 1365000.0
    assert round(value['not_for'], 4) == 1000000.0
    assert value['ccy_dom'] == 'USD'
    assert value['ccy_for'] == 'EUR'

    delta = call_option.delta(
        valuation_date,
        spot_fx_rate,
        dom_discount_curve,
        for_discount_curve,
        model)

    assert round(delta['pips_spot_delta'], 4) == 0.3671
    assert round(delta['pips_fwd_delta'], 4) == 0.3859
    assert round(delta['pct_spot_delta_prem_adj'], 4) == 0.3208
    assert round(delta['pct_fwd_delta_prem_adj'], 4) == 0.3373


def test_value_mc():
    #   Example from Hull 4th edition page 284
    valuation_date = Date(1, 1, 2015)
    expiry_date = valuation_date.add_months(4)
    spot_fx_rate = 1.60
    volatility = 0.1411
    dom_interest_rate = 0.08
    forInterestRate = 0.11
    model = BlackScholes(volatility)
    dom_discount_curve = DiscountCurveFlat(valuation_date, dom_interest_rate)
    for_discount_curve = DiscountCurveFlat(valuation_date, forInterestRate)
    num_paths = 100000

    strike_fx_rate = 1.6

    call_option = FXVanillaOption(expiry_date,
                                  strike_fx_rate,
                                  "EURUSD",
                                  OptionTypes.EUROPEAN_CALL,
                                  1000000,
                                  "USD")

    value_mc = call_option.value_mc(
        valuation_date,
        spot_fx_rate,
        dom_discount_curve,
        for_discount_curve,
        model,
        num_paths)

    assert round(value_mc, 4) == 0.0429

    put_option = FXVanillaOption(expiry_date,
                                 strike_fx_rate,
                                 "EURUSD",
                                 OptionTypes.EUROPEAN_PUT,
                                 1000000,
                                 "USD")

    value_mc = put_option.value_mc(
        valuation_date,
        spot_fx_rate,
        dom_discount_curve,
        for_discount_curve,
        model,
        num_paths)

    assert round(value_mc, 4) == 0.0582


def test_vega_theta():
    #   Example from Hull 4th edition page 284
    valuation_date = Date(1, 1, 2015)
    expiry_date = valuation_date.add_months(4)
    spot_fx_rate = 1.60
    volatility = 0.1411
    dom_interest_rate = 0.08
    forInterestRate = 0.11
    model = BlackScholes(volatility)
    dom_discount_curve = DiscountCurveFlat(valuation_date, dom_interest_rate)
    for_discount_curve = DiscountCurveFlat(valuation_date, forInterestRate)

    strike_fx_rate = 1.6

    call_option = FXVanillaOption(expiry_date,
                                  strike_fx_rate,
                                  "EURUSD",
                                  OptionTypes.EUROPEAN_CALL,
                                  1000000,
                                  "USD")

    vega = call_option.vega(
        valuation_date,
        spot_fx_rate,
        dom_discount_curve,
        for_discount_curve,
        model)

    assert round(vega, 4) == 0.3518

    theta = call_option.theta(
        valuation_date,
        spot_fx_rate,
        dom_discount_curve,
        for_discount_curve,
        model)

    assert round(theta, 4) == -0.0504
