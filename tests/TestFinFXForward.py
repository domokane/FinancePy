# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.fx.FinFXForward import FinFXForward
from financepy.market.curves.FinFlatCurve import FinFlatCurve

from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinCalendar import FinCalendarTypes
from financepy.market.curves.FinLiborCurve import FinLiborCurve
from financepy.products.libor.FinLiborDeposit import FinLiborDeposit

from financepy.finutils.FinDate import FinDate
import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinFXForward():

    #  https://stackoverflow.com/questions/48778712
    #  /fx-vanilla-call-price-in-quantlib-doesnt-match-bloomberg

    valueDate = FinDate(13, 2, 2018)
    expiryDate = valueDate.addMonths(12)
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
    settlementDate = valueDate.addWorkDays(spotDays)
    maturityDate = settlementDate.addMonths(12)
    notional = 100.0
    calendarType = FinCalendarTypes.TARGET

    depos = []
    fras = []
    swaps = []
    depositRate = ccy1InterestRate
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate,
                           FinDayCountTypes.ACT_360, notional, calendarType)
    depos.append(depo)
    forDiscountCurve = FinLiborCurve(forName, settlementDate, depos, fras, swaps)

    depos = []
    fras = []
    swaps = []
    depositRate = ccy2InterestRate
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate,
                           FinDayCountTypes.ACT_360, notional, calendarType)
    depos.append(depo)
    domDiscountCurve = FinLiborCurve(domName, settlementDate, depos, fras, swaps)

    notional = 100.0
    notionalCurrency = forName

    fxForward = FinFXForward(expiryDate,
                             strikeFXRate,
                             currencyPair,
                             notional,
                             notionalCurrency)

    testCases.header("SPOT FX", "FX FWD", "VALUE_BS")

    fwdValue = fxForward.value(valueDate, spotFXRate,
                               domDiscountCurve, forDiscountCurve)

    fwdFXRate = fxForward.forward(valueDate, spotFXRate,
                                  domDiscountCurve,
                                  forDiscountCurve)

    testCases.print(spotFXRate, fwdFXRate, fwdValue)

###############################################################################


test_FinFXForward()
testCases.compareTestCases()
