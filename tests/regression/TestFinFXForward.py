##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import add_fp_to_path

from financepy.utils.date import Date
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.market.curves.ibor_single_curve import IborSingleCurve
from financepy.utils.calendar import CalendarTypes
from financepy.utils.day_count import DayCountTypes

from financepy.products.fx.fx_forward import FXForward, fx_forward_rate
from financepy.products.fx.fx_rate import FXRate

from FinTestCases import FinTestCases, global_test_case_mode


test_cases = FinTestCases(__file__, global_test_case_mode)


##############################################################################


def test_fin_fx_forward():

    # Valuation date and currency pair.
    value_dt = Date(13, 2, 2018)
    currency_pair = "EURUSD"

    # FX rates are quoted USD per EUR.
    spot_fx_rate = FXRate(currency_pair, 1.300)
    strike_fx_rate = FXRate(currency_pair, 1.365)

    # EUR and USD interest rates.
    eur_interest_rate = 0.02
    usd_interest_rate = 0.05

    # Spot and forward settlement dates.
    spot_dt = spot_fx_rate.spot_date(value_dt)
    settlement_dt = spot_dt.add_months(12)

    notional = 100.0
    notional_currency = "EUR"

    ##########################################################################
    # EUR discount curve: first currency in EURUSD.
    ##########################################################################

    euro_cal_type = CalendarTypes.TARGET

    eur_depo = IborDeposit(
        value_dt,
        settlement_dt,
        eur_interest_rate,
        DayCountTypes.ACT_360,
        notional,
        euro_cal_type,
    )

    foreign_curve = IborSingleCurve(
        value_dt,
        [eur_depo],
        [],
        [],
    )

    ##########################################################################
    # USD discount curve: second currency in EURUSD.
    ##########################################################################

    usd_cal_type = CalendarTypes.NEW_YORK

    usd_depo = IborDeposit(
        value_dt,
        settlement_dt,
        usd_interest_rate,
        DayCountTypes.ACT_360,
        notional,
        usd_cal_type,
    )

    domestic_curve = IborSingleCurve(
        value_dt,
        [usd_depo],
        [],
        [],
    )

    ##########################################################################
    # Buy EUR, sell USD.
    ##########################################################################

    fx_fwd_buy = FXForward(
        settlement_dt=settlement_dt,
        strike_fx_rate=strike_fx_rate,
        notional=notional,
        notional_ccy=notional_currency,
        buy_ccy="EUR",
    )

    buy_value = fx_fwd_buy.value(
        value_dt,
        spot_fx_rate,
        domestic_curve,
        foreign_curve,
    )

    fair_fwd_rate = fx_forward_rate(
        value_dt,
        settlement_dt,
        spot_fx_rate,
        domestic_curve,
        foreign_curve,
    )

    test_cases.header(
        "SPOT FX",
        "FX FORWARD",
        "VALUE",
    )

    test_cases.print(
        spot_fx_rate,
        fair_fwd_rate,
        buy_value,
    )

    ##########################################################################
    # Sell EUR, buy USD.
    ##########################################################################

    fx_fwd_sell = FXForward(
        settlement_dt=settlement_dt,
        strike_fx_rate=strike_fx_rate,
        notional=notional,
        notional_ccy=notional_currency,
        buy_ccy="USD",
    )

    sell_value = fx_fwd_sell.value(
        value_dt,
        spot_fx_rate,
        domestic_curve,
        foreign_curve,
    )

    test_cases.header(
        "BUY VALUE",
        "SELL VALUE",
        "SUM",
    )

    test_cases.print(
        buy_value["USD"],
        sell_value["USD"],
        buy_value["USD"] + sell_value["USD"],
    )

    assert abs(
        buy_value["USD"] + sell_value["USD"]
    ) < 1e-10

    assert abs(
        buy_value["EUR"] + sell_value["EUR"]
    ) < 1e-10

    ##########################################################################
    # Forward struck at the fair forward rate.
    # Its initial value should be zero.
    ##########################################################################

    fx_fwd_fair = FXForward(
        settlement_dt=settlement_dt,
        strike_fx_rate=fair_fwd_rate,
        notional=notional,
        notional_ccy="EUR",
        buy_ccy="EUR",
    )

    fair_value = fx_fwd_fair.value(
        value_dt,
        spot_fx_rate,
        domestic_curve,
        foreign_curve,
    )

    test_cases.header(
        "FAIR FORWARD RATE",
        "NPV CCY1",
        "NPV CCY2",
    )

    test_cases.print(
        fair_fwd_rate,
        fair_value["EUR"],
        fair_value["USD"],
    )

    assert abs(fair_value["EUR"]) < 1e-10
    assert abs(fair_value["USD"]) < 1e-10

    ##########################################################################
    # Spot settlement test.
    # A new transaction struck at the market spot rate has zero value.
    ##########################################################################

    fx_spot_equivalent = FXForward(
        settlement_dt=spot_dt,
        strike_fx_rate=spot_fx_rate,
        notional=notional,
        notional_ccy="EUR",
        buy_ccy="EUR",
    )

    spot_value = fx_spot_equivalent.value(
        value_dt,
        spot_fx_rate,
        domestic_curve,
        foreign_curve,
    )

    spot_forward_rate = fx_forward_rate(
        value_dt,
        spot_dt,
        spot_fx_rate,
        domestic_curve,
        foreign_curve,
    )

    test_cases.header(
        "MARKET SPOT",
        "SPOT FORWARD",
        "SPOT NPV",
    )

    test_cases.print(
        spot_fx_rate,
        spot_forward_rate,
        spot_value["USD"],
    )

    assert abs(
        spot_forward_rate.rate - spot_fx_rate.rate
    ) < 1e-12

    assert abs(spot_value["USD"]) < 1e-10

    ##########################################################################
    # Contractual cash-flow checks.
    ##########################################################################

    assert fx_fwd_buy.cashflow_for == 100.0
    assert abs(fx_fwd_buy.cashflow_dom + 136.5) < 1e-12

    assert fx_fwd_sell.cashflow_for == -100.0
    assert abs(fx_fwd_sell.cashflow_dom - 136.5) < 1e-12

##############################################################################


test_fin_fx_forward()
test_cases.compare_test_cases()
