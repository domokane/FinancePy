###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.date import Date
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.utils.calendar import CalendarTypes
from financepy.utils.day_count import DayCountTypes
from financepy.products.fx.fx_forward import FXForward


def test_FinFXForward():
    #  https://stackoverflow.com/questions/48778712
    #  /fx-vanilla-call-price-in-quantlib-doesnt-match-bloomberg

    value_dt = Date(13, 2, 2018)
    expiry_dt = value_dt.add_months(12)
    # Forward is on EURUSD which is expressed as number of USD per EUR
    # ccy1 = EUR and ccy2 = USD
    for_name = "EUR"
    dom_name = "USD"
    currency_pair = for_name + dom_name  # Always ccy1ccy2
    spot_fx_rate = 1.300  # USD per EUR
    strike_fx_rate = 1.365  # USD per EUR
    ccy1InterestRate = 0.02  # USD Rates
    ccy2InterestRate = 0.05  # EUR rates

    spot_days = 0
    settle_dt = value_dt.add_weekdays(spot_days)
    maturity_dt = settle_dt.add_months(12)
    notional = 100.0
    cal_type = CalendarTypes.TARGET

    depos = []
    fras = []
    swaps = []
    deposit_rate = ccy1InterestRate
    depo = IborDeposit(
        settle_dt,
        maturity_dt,
        deposit_rate,
        DayCountTypes.ACT_360,
        notional,
        cal_type,
    )
    depos.append(depo)
    foreign_curve = IborSingleCurve(value_dt, depos, fras, swaps)

    depos = []
    fras = []
    swaps = []
    deposit_rate = ccy2InterestRate
    depo = IborDeposit(
        settle_dt,
        maturity_dt,
        deposit_rate,
        DayCountTypes.ACT_360,
        notional,
        cal_type,
    )
    depos.append(depo)
    domestic_curve = IborSingleCurve(value_dt, depos, fras, swaps)

    notional = 100.0
    notional_currency = for_name

    fx_fwd = FXForward(
        expiry_dt, strike_fx_rate, currency_pair, notional, notional_currency
    )

    fwd_value = fx_fwd.value(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve
    )

    fwd_fx_rate = fx_fwd.forward(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve
    )

    assert round(fwd_fx_rate, 4) == 1.3388

    assert round(fwd_value["value"], 4) == -2.4978
    assert round(fwd_value["cash_dom"], 4) == -249.7797
    assert round(fwd_value["cash_for"], 4) == -192.1382
    assert fwd_value["not_dom"] == 136.5
    assert fwd_value["not_for"] == 100.0
    assert fwd_value["ccy_dom"] == "USD"
    assert fwd_value["ccy_for"] == "EUR"
