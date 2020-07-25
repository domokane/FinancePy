# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 16:23:12 2019

@author: Dominic
"""

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.products.libor.FinLiborCurve import FinLiborCurve
from financepy.finutils.FinOptionTypes import FinOptionExerciseTypes

from financepy.products.libor.FinLiborDeposit import FinLiborDeposit
from financepy.products.libor.FinLiborSwap import FinLiborSwap
from financepy.products.libor.FinLiborSwaption import FinLiborSwaptionTypes
from financepy.products.libor.FinLiborSwaption import FinLiborSwaption
from financepy.products.libor.FinLiborBermudanSwaption import FinLiborBermudanSwaption
from financepy.models.FinModelBlack import FinModelBlack
from financepy.models.FinModelRatesBK import FinModelRatesBK
from financepy.models.FinModelRatesHW import FinModelRatesHW
from financepy.models.FinModelRatesBDT import FinModelRatesBDT
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat

testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinLiborBermudanSwaptionBKModel():
    ''' Replicate examples in paper by Leif Andersen which can be found at
    file:///C:/Users/Dominic/Downloads/SSRN-id155208.pdf '''

    import time

    valuationDate = FinDate(2011, 1, 1)
    settlementDate = valuationDate
    exerciseDate = settlementDate.addYears(1)
    swapMaturityDate = settlementDate.addYears(4)

    swapFixedCoupon = 0.060
    swapFixedFrequencyType = FinFrequencyTypes.SEMI_ANNUAL
    swapFixedDayCountType = FinDayCountTypes.ACT_365_ISDA

    liborCurve = FinDiscountCurveFlat(valuationDate,
                                      0.06,
                                      FinFrequencyTypes.SEMI_ANNUAL)

    start = time.time()

    model = FinModelBlack(0.25)

    swaptionType = FinLiborSwaptionTypes.PAYER
    europeanSwaptionPayer = FinLiborSwaption(settlementDate,
                                             exerciseDate,
                                             swapMaturityDate,
                                             swaptionType,
                                             swapFixedCoupon,
                                             swapFixedFrequencyType,
                                             swapFixedDayCountType)

    value = europeanSwaptionPayer.value(settlementDate, liborCurve, model)
    print("European Black Payer Value:", value)

    swaptionType = FinLiborSwaptionTypes.RECEIVER
    europeanSwaptionReceiver = FinLiborSwaption(settlementDate,
                                                exerciseDate,
                                                swapMaturityDate,
                                                swaptionType,
                                                swapFixedCoupon,
                                                swapFixedFrequencyType,
                                                swapFixedDayCountType)

    value = europeanSwaptionReceiver.value(valuationDate, liborCurve, model)
    print("Euroopean Black Recev Value:", value)

    swaptionType = FinLiborSwaptionTypes.PAYER
    exerciseType = FinOptionExerciseTypes.BERMUDAN

    # Andersen used BDT with constant short-rate volatility
    sigma = 0.2012
    numTimeSteps = 5
    model = FinModelRatesBDT(sigma, numTimeSteps)

    bermudanSwaptionPayer = FinLiborBermudanSwaption(settlementDate,
                                                     exerciseDate,
                                                     swapMaturityDate,
                                                     swaptionType,
                                                     exerciseType,
                                                     swapFixedCoupon,
                                                     swapFixedFrequencyType,
                                                     swapFixedDayCountType)

    value = bermudanSwaptionPayer.value(valuationDate, liborCurve, model)

    print("Bermudan Payer Value:", value)

    end = time.time()

    testCases.header("TIME")
    testCases.print(end - start)

##########################################################################


test_FinLiborBermudanSwaptionBKModel()

testCases.compareTestCases()
