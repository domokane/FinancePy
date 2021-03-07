###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.products.fx.fx_forward import FinFXForward
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinCalendar import FinCalendarTypes
from financepy.products.rates.FinIborSingleCurve import FinIborSingleCurve
from financepy.products.rates.FinIborDeposit import FinIborDeposit
from financepy.finutils.FinDate import FinDate

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinFXForward():

    #  https://stackoverflow.com/questions/48778712
    #  /fx-vanilla-call-price-in-quantlib-doesnt-match-bloomberg

    valuationDate = FinDate(13, 2, 2018)
    expiryDate = valuationDate.addMonths(12)
    # Forward is on EURUSD which is expressed as number of USD per EUR
    # ccy1 = EUR and ccy2 = USD
    forName = "EUR"
    domName = "USD"
    currency_pair = forName + domName  # Always ccy1ccy2
    spot_fx_rate = 1.300  # USD per EUR
    strike_fx_rate = 1.365  # USD per EUR
    ccy1InterestRate = 0.02  # USD Rates
    ccy2InterestRate = 0.05  # EUR rates

    ###########################################################################

    spot_days = 0
    settlementDate = valuationDate.addWeekDays(spot_days)
    maturityDate = settlementDate.addMonths(12)
    notional = 100.0
    calendarType = FinCalendarTypes.TARGET

    depos = []
    fras = []
    swaps = []
    depositRate = ccy1InterestRate
    depo = FinIborDeposit(settlementDate, maturityDate, depositRate,
                           FinDayCountTypes.ACT_360, notional, calendarType)
    depos.append(depo)
    for_discount_curve = FinIborSingleCurve(valuationDate, depos, fras, swaps)

    depos = []
    fras = []
    swaps = []
    depositRate = ccy2InterestRate
    depo = FinIborDeposit(settlementDate, maturityDate, depositRate,
                           FinDayCountTypes.ACT_360, notional, calendarType)
    depos.append(depo)
    dom_discount_curve = FinIborSingleCurve(valuationDate, depos, fras, swaps)

    notional = 100.0
    notional_currency = forName

    fxForward = FinFXForward(expiryDate,
                             strike_fx_rate,
                             currency_pair,
                             notional,
                             notional_currency)

    testCases.header("SPOT FX", "FX FWD", "VALUE_BS")

    fwdValue = fxForward.value(valuationDate, spot_fx_rate,
                               dom_discount_curve, for_discount_curve)

    fwdFXRate = fxForward.forward(valuationDate, spot_fx_rate,
                                  dom_discount_curve,
                                  for_discount_curve)

    testCases.print(spot_fx_rate, fwdFXRate, fwdValue)

###############################################################################


test_FinFXForward()
testCases.compareTestCases()
