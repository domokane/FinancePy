###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinOptionTypes import FinOptionExerciseTypes
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

    valuationDate = FinDate(1, 1, 2011)
    settlementDate = valuationDate
    exerciseDate = settlementDate.addYears(1)
    swapMaturityDate = settlementDate.addYears(4)

    swapFixedCoupon = 0.060
    swapFixedFrequencyType = FinFrequencyTypes.SEMI_ANNUAL
    swapFixedDayCountType = FinDayCountTypes.ACT_365F

    liborCurve = FinDiscountCurveFlat(valuationDate,
                                      0.06,
                                      FinFrequencyTypes.SEMI_ANNUAL)

    ###########################################################################
    # BLACK'S MODEL
    ###########################################################################

    model = FinModelBlack(0.25)

    swaptionType = FinLiborSwaptionTypes.PAYER
    europeanSwaptionPay = FinLiborSwaption(settlementDate,
                                           exerciseDate,
                                           swapMaturityDate,
                                           swaptionType,
                                           swapFixedCoupon,
                                           swapFixedFrequencyType,
                                           swapFixedDayCountType)

    value = europeanSwaptionPay.value(settlementDate, liborCurve, model)
    print("EUROPEAN BLACK PAY VALUE:", value)

    swaptionType = FinLiborSwaptionTypes.RECEIVER
    europeanSwaptionRec = FinLiborSwaption(settlementDate,
                                           exerciseDate,
                                           swapMaturityDate,
                                           swaptionType,
                                           swapFixedCoupon,
                                           swapFixedFrequencyType,
                                           swapFixedDayCountType)

    value = europeanSwaptionRec.value(settlementDate, liborCurve, model)
    print("EUROPEAN BLACK REC VALUE:", value)

    ###########################################################################
    # BK MODEL
    ###########################################################################

    # Andersen used BDT with constant short-rate volatility
    sigma = 0.2012
    a = 0.01
    numTimeSteps = 100
    model = FinModelRatesBK(sigma, a, numTimeSteps)

    value = europeanSwaptionPay.value(valuationDate, liborCurve, model)
    print("EUROPEAN BK PAY VALUE:", value)

    value = europeanSwaptionRec.value(valuationDate, liborCurve, model)
    print("EUROPEAN BK REC VALUE:", value)

    ###########################################################################
    # BDT MODEL
    ###########################################################################

    # Andersen used BDT with constant short-rate volatility
    sigma = 0.2012
    numTimeSteps = 100
    model = FinModelRatesBDT(sigma, numTimeSteps)

    value = europeanSwaptionPay.value(valuationDate, liborCurve, model)
    print("EUROPEAN BK PAY VALUE:", value)

    value = europeanSwaptionRec.value(valuationDate, liborCurve, model)
    print("EUROPEAN BK REC VALUE:", value)

    ###########################################################################

    swaptionType = FinLiborSwaptionTypes.PAYER
    exerciseType = FinOptionExerciseTypes.BERMUDAN

    bermudanSwaptionPay = FinLiborBermudanSwaption(settlementDate,
                                                   exerciseDate,
                                                   swapMaturityDate,
                                                   swaptionType,
                                                   exerciseType,
                                                   swapFixedCoupon,
                                                   swapFixedFrequencyType,
                                                   swapFixedDayCountType)

    swaptionType = FinLiborSwaptionTypes.RECEIVER
    exerciseType = FinOptionExerciseTypes.BERMUDAN

    bermudanSwaptionRec = FinLiborBermudanSwaption(settlementDate,
                                                   exerciseDate,
                                                   swapMaturityDate,
                                                   swaptionType,
                                                   exerciseType,
                                                   swapFixedCoupon,
                                                   swapFixedFrequencyType,
                                                   swapFixedDayCountType)

    ###########################################################################
    # USING BK MODEL
    ###########################################################################

    sigma = 0.2012
    a = 0.1
    numTimeSteps = 100
    model = FinModelRatesBK(sigma, a, numTimeSteps)

    value = bermudanSwaptionPay.value(valuationDate, liborCurve, model)
    print("Bermudan BK PAY Value:", value)
    value = bermudanSwaptionRec.value(valuationDate, liborCurve, model)
    print("Bermudan BK REC Value:", value)

    ###########################################################################
    # BDT MODEL
    ###########################################################################

    # Andersen used BDT with constant short-rate volatility
    sigma = 0.2012
    numTimeSteps = 100
    model = FinModelRatesBDT(sigma, numTimeSteps)

    value = bermudanSwaptionPay.value(valuationDate, liborCurve, model)
    print("Bermudan BDT PAY Value:", value)
    value = bermudanSwaptionRec.value(valuationDate, liborCurve, model)
    print("Bermudan BDT REC Value:", value)

    ###########################################################################
    # USING HW MODEL
    ###########################################################################

    sigma = 0.01
    a = 0.1
    numTimeSteps = 100
    model = FinModelRatesHW(sigma, a, numTimeSteps)

    value = bermudanSwaptionPay.value(valuationDate, liborCurve, model)
    print("Bermudan HW PAY Value:", value)
    value = bermudanSwaptionRec.value(valuationDate, liborCurve, model)
    print("Bermudan HW REC Value:", value)

##########################################################################


test_FinLiborBermudanSwaptionBKModel()

testCases.compareTestCases()
