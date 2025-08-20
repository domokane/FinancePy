########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.date import Date
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.utils.calendar import CalendarTypes
from financepy.utils.day_count import DayCountTypes
from financepy.products.fx.fx_forward import FXForward
import sys

sys.path.append("..")


test_cases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


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
    ccy1_interest_rate = 0.02  # USD Rates
    ccy2_interest_rate = 0.05  # EUR rates

    ###########################################################################

    spot_days = 0
    settle_dt = value_dt.add_weekdays(spot_days)
    maturity_dt = settle_dt.add_months(12)
    notional = 100.0
    cal_type = CalendarTypes.TARGET

    depos = []
    fras = []
    swaps = []
    deposit_rate = ccy1_interest_rate
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
    deposit_rate = ccy2_interest_rate
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

    test_cases.header("SPOT FX", "FX FWD", "VALUE_BS")

    fwd_value = fx_fwd.value(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve
    )

    fwd_fx_rate = fx_fwd.forward(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve
    )

    test_cases.print(spot_fx_rate, fwd_fx_rate, fwd_value)


########################################################################################


test_FinFXForward()
test_cases.compareTestCases()
