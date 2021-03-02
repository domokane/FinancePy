###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.products.fx.FinFXForward import FinFXForward
from financepy.utils.day_count import DayCountTypes
from financepy.utils.calendar import CalendarTypes
from financepy.products.rates.FinIborSingleCurve import IborSingleCurve
from financepy.products.rates.FinIborDeposit import FinIborDeposit
from financepy.utils.date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinFXForward():

    #  https://stackoverflow.com/questions/48778712
    #  /fx-vanilla-call-price-in-quantlib-doesnt-match-bloomberg

    valuation_date = Date(13, 2, 2018)
    expiry_date = valuation_date.addMonths(12)
    # Forward is on EURUSD which is expressed as number of USD per EUR
    # ccy1 = EUR and ccy2 = USD
    forName = "EUR"
    domName = "USD"
    currencyPair = forName + domName  # Always ccy1ccy2
    spotFXRate = 1.300  # USD per EUR
    strikeFXRate = 1.365  # USD per EUR
    ccy1InterestRate = 0.02  # USD Rates
    ccy2InterestRate = 0.05  # EUR rates

    ###########################################################################

    spotDays = 0
    settlement_date = valuation_date.addWeekDays(spotDays)
    maturity_date = settlement_date.addMonths(12)
    notional = 100.0
    calendar_type = CalendarTypes.TARGET

    depos = []
    fras = []
    swaps = []
    deposit_rate = ccy1InterestRate
    depo = FinIborDeposit(settlement_date, maturity_date, deposit_rate,
                           DayCountTypes.ACT_360, notional, calendar_type)
    depos.append(depo)
    forDiscountCurve = IborSingleCurve(valuation_date, depos, fras, swaps)

    depos = []
    fras = []
    swaps = []
    deposit_rate = ccy2InterestRate
    depo = FinIborDeposit(settlement_date, maturity_date, deposit_rate,
                           DayCountTypes.ACT_360, notional, calendar_type)
    depos.append(depo)
    domDiscountCurve = IborSingleCurve(valuation_date, depos, fras, swaps)

    notional = 100.0
    notionalCurrency = forName

    fxForward = FinFXForward(expiry_date,
                             strikeFXRate,
                             currencyPair,
                             notional,
                             notionalCurrency)

    testCases.header("SPOT FX", "FX FWD", "VALUE_BS")

    fwdValue = fxForward.value(valuation_date, spotFXRate,
                               domDiscountCurve, forDiscountCurve)

    fwdFXRate = fxForward.forward(valuation_date, spotFXRate,
                                  domDiscountCurve,
                                  forDiscountCurve)

    testCases.print(spotFXRate, fwdFXRate, fwdValue)

###############################################################################


test_FinFXForward()
testCases.compareTestCases()
