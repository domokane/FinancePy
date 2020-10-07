###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..//financepy")

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinGlobalTypes import FinSwapTypes
from financepy.finutils.FinGlobalTypes import FinExerciseTypes
from financepy.products.libor.FinLiborSwaption import FinLiborSwaption
from financepy.products.libor.FinLiborSwap import FinLiborSwap

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
                                      0.0625,
                                      FinFrequencyTypes.SEMI_ANNUAL)

    fwdPayerSwap = FinLiborSwap(exerciseDate,
                                swapMaturityDate,
                                FinSwapTypes.PAYER,
                                swapFixedCoupon,
                                swapFixedFrequencyType,
                                swapFixedDayCountType)

    fwdSwapValue = fwdPayerSwap.value(settlementDate, liborCurve, liborCurve)

    testCases.header("LABEL", "VALUE")
    testCases.print("FWD SWAP VALUE", fwdSwapValue)

    # fwdPayerSwap.printFixedLegPV()

    # Now we create the European swaptions
    swapType = FinSwapTypes.PAYER
    europeanSwaptionPay = FinLiborSwaption(settlementDate,
                                           exerciseDate,
                                           swapMaturityDate,
                                           swapType,
                                           swapFixedCoupon,
                                           swapFixedFrequencyType,
                                           swapFixedDayCountType)

    swapType = FinSwapTypes.RECEIVER
    europeanSwaptionRec = FinLiborSwaption(settlementDate,
                                           exerciseDate,
                                           swapMaturityDate,
                                           swapType,
                                           swapFixedCoupon,
                                           swapFixedFrequencyType,
                                           swapFixedDayCountType)
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    # BLACK'S MODEL
    ###########################################################################
    ###########################################################################
    ###########################################################################

    testCases.banner("======= ZERO VOLATILITY ========")
    model = FinModelBlack(0.0000001)
    testCases.print("Black Model", model._volatility)

    valuePay = europeanSwaptionPay.value(settlementDate, liborCurve, model)
    testCases.print("EUROPEAN BLACK PAY VALUE ZERO VOL:", valuePay)

    valueRec = europeanSwaptionRec.value(settlementDate, liborCurve, model)
    testCases.print("EUROPEAN BLACK REC VALUE ZERO VOL:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAYER MINUS RECEIVER :", payRec)

    testCases.banner("======= 20%% BLACK VOLATILITY ========")

    model = FinModelBlack(0.20)
    testCases.print("Black Model", model._volatility)

    valuePay = europeanSwaptionPay.value(settlementDate, liborCurve, model)
    testCases.print("EUROPEAN BLACK PAY VALUE:", valuePay)

    valueRec = europeanSwaptionRec.value(settlementDate, liborCurve, model)
    testCases.print("EUROPEAN BLACK REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAYER MINUS RECEIVER :", payRec)

    ###########################################################################
    ###########################################################################
    ###########################################################################
    # BK MODEL
    ###########################################################################
    ###########################################################################
    ###########################################################################

    testCases.banner("=======================================================")
    testCases.banner("=======================================================")
    testCases.banner("==================== BK MODEL =========================")
    testCases.banner("=======================================================")
    testCases.banner("=======================================================")

    testCases.banner("======= 0% VOLATILITY EUROPEAN SWAPTION BK MODEL ======")

    # Used BK with constant short-rate volatility
    sigma = 0.00001
    a = 0.01
    numTimeSteps = 200
    model = FinModelRatesBK(sigma, a, numTimeSteps)

    valuePay = europeanSwaptionPay.value(valuationDate, liborCurve, model)
    testCases.print("EUROPEAN BK PAY VALUE:", valuePay)

    valueRec = europeanSwaptionRec.value(valuationDate, liborCurve, model)
    testCases.print("EUROPEAN BK REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAYER MINUS RECEIVER :", payRec)


    testCases.banner("======= 20% VOLATILITY EUROPEAN SWAPTION BK MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.20
    a = 0.01
    model = FinModelRatesBK(sigma, a, numTimeSteps)

    testCases.banner("BK MODEL SWAPTION CLASS EUROPEAN EXERCISE")

    valuePay = europeanSwaptionPay.value(valuationDate, liborCurve, model)
    testCases.print("EUROPEAN BK PAY VALUE:", valuePay)

    valueRec = europeanSwaptionRec.value(valuationDate, liborCurve, model)
    testCases.print("EUROPEAN BK REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAYER MINUS RECEIVER :", payRec)
    
    ###########################################################################

    # Now we create the Bermudan swaptions but only allow European exercise
    swapType = FinSwapTypes.PAYER
    exerciseType = FinExerciseTypes.EUROPEAN

    bermudanSwaptionPay = FinLiborBermudanSwaption(settlementDate,
                                                   exerciseDate,
                                                   swapMaturityDate,
                                                   swapType,
                                                   exerciseType,
                                                   swapFixedCoupon,
                                                   swapFixedFrequencyType,
                                                   swapFixedDayCountType)

    swapType = FinSwapTypes.RECEIVER
    exerciseType = FinExerciseTypes.EUROPEAN

    bermudanSwaptionRec = FinLiborBermudanSwaption(settlementDate,
                                                   exerciseDate,
                                                   swapMaturityDate,
                                                   swapType,
                                                   exerciseType,
                                                   swapFixedCoupon,
                                                   swapFixedFrequencyType,
                                                   swapFixedDayCountType)
   
    testCases.banner("======= 0% VOLATILITY BERMUDAN SWAPTION EUROPEAN EXERCISE BK MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.000001
    a = 0.01
    model = FinModelRatesBK(sigma, a, numTimeSteps)

    testCases.banner("BK MODEL BERMUDAN SWAPTION CLASS EUROPEAN EXERCISE")
    valuePay = bermudanSwaptionPay.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN BK PAY VALUE:", valuePay)

    valueRec = bermudanSwaptionRec.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN BK REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAYER MINUS RECEIVER :", payRec)

    testCases.banner("======= 20% VOLATILITY BERMUDAN SWAPTION EUROPEAN EXERCISE BK MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.2
    a = 0.01
    model = FinModelRatesBK(sigma, a, numTimeSteps)

    testCases.banner("BK MODEL BERMUDAN SWAPTION CLASS EUROPEAN EXERCISE")
    valuePay = bermudanSwaptionPay.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN BK PAY VALUE:", valuePay)

    valueRec = bermudanSwaptionRec.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN BK REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAYER MINUS RECEIVER :", payRec)
    
    ###########################################################################
    # Now we create the Bermudan swaptions but allow Bermudan exercise
    ###########################################################################

    swapType = FinSwapTypes.PAYER
    exerciseType = FinExerciseTypes.BERMUDAN

    bermudanSwaptionPay = FinLiborBermudanSwaption(settlementDate,
                                                   exerciseDate,
                                                   swapMaturityDate,
                                                   swapType,
                                                   exerciseType,
                                                   swapFixedCoupon,
                                                   swapFixedFrequencyType,
                                                   swapFixedDayCountType)

    swapType = FinSwapTypes.RECEIVER
    exerciseType = FinExerciseTypes.BERMUDAN

    bermudanSwaptionRec = FinLiborBermudanSwaption(settlementDate,
                                                   exerciseDate,
                                                   swapMaturityDate,
                                                   swapType,
                                                   exerciseType,
                                                   swapFixedCoupon,
                                                   swapFixedFrequencyType,
                                                   swapFixedDayCountType)

    testCases.banner("======= ZERO VOLATILITY BERMUDAN SWAPTION BERMUDAN EXERCISE BK MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.000001
    a = 0.01
    model = FinModelRatesBK(sigma, a, numTimeSteps)

    testCases.banner("BK MODEL BERMUDAN SWAPTION CLASS BERMUDAN EXERCISE")
    valuePay = bermudanSwaptionPay.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN BK PAY VALUE:", valuePay)

    valueRec = bermudanSwaptionRec.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN BK REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAYER MINUS RECEIVER :", payRec)

    testCases.banner("======= 20% VOLATILITY BERMUDAN SWAPTION BERMUDAN EXERCISE BK MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.20
    a = 0.01
    model = FinModelRatesBK(sigma, a, numTimeSteps)

    testCases.banner("BK MODEL BERMUDAN SWAPTION CLASS BERMUDAN EXERCISE")
    valuePay = bermudanSwaptionPay.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN BK PAY VALUE:", valuePay)

    valueRec = bermudanSwaptionRec.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN BK REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAYER MINUS RECEIVER :", payRec)
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    # BDT MODEL
    ###########################################################################
    ###########################################################################
    ###########################################################################

    testCases.banner("=======================================================")
    testCases.banner("=======================================================")
    testCases.banner("======================= BDT MODEL =====================")
    testCases.banner("=======================================================")
    testCases.banner("=======================================================")

    testCases.banner("====== 0% VOLATILITY EUROPEAN SWAPTION BDT MODEL ======")

    # Used BK with constant short-rate volatility
    sigma = 0.00001
    numTimeSteps = 200
    model = FinModelRatesBDT(sigma, numTimeSteps)

    valuePay = europeanSwaptionPay.value(valuationDate, liborCurve, model)
    testCases.print("EUROPEAN BDT PAY VALUE:", valuePay)

    valueRec = europeanSwaptionRec.value(valuationDate, liborCurve, model)
    testCases.print("EUROPEAN BDT REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAYER MINUS RECEIVER :", payRec)

    testCases.banner("===== 20% VOLATILITY EUROPEAN SWAPTION BDT MODEL ======")

    # Used BK with constant short-rate volatility
    sigma = 0.20
    a = 0.01
    model = FinModelRatesBDT(sigma, numTimeSteps)

    testCases.banner("BDT MODEL SWAPTION CLASS EUROPEAN EXERCISE")

    valuePay = europeanSwaptionPay.value(valuationDate, liborCurve, model)
    testCases.print("EUROPEAN BDT PAY VALUE:", valuePay)

    valueRec = europeanSwaptionRec.value(valuationDate, liborCurve, model)
    testCases.print("EUROPEAN BDT REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAYER MINUS RECEIVER :", payRec)

    ###########################################################################

    # Now we create the Bermudan swaptions but only allow European exercise
    swapType = FinSwapTypes.PAYER
    exerciseType = FinExerciseTypes.EUROPEAN

    bermudanSwaptionPay = FinLiborBermudanSwaption(settlementDate,
                                                   exerciseDate,
                                                   swapMaturityDate,
                                                   swapType,
                                                   exerciseType,
                                                   swapFixedCoupon,
                                                   swapFixedFrequencyType,
                                                   swapFixedDayCountType)

    swapType = FinSwapTypes.RECEIVER
    bermudanSwaptionRec = FinLiborBermudanSwaption(settlementDate,
                                                   exerciseDate,
                                                   swapMaturityDate,
                                                   swapType,
                                                   exerciseType,
                                                   swapFixedCoupon,
                                                   swapFixedFrequencyType,
                                                   swapFixedDayCountType)
   
    testCases.banner("======= 0% VOLATILITY BERMUDAN SWAPTION EUROPEAN EXERCISE BDT MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.000001
    model = FinModelRatesBDT(sigma, numTimeSteps)

    testCases.banner("BK MODEL BERMUDAN SWAPTION CLASS EUROPEAN EXERCISE")
    valuePay = bermudanSwaptionPay.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN BDT PAY VALUE:", valuePay)

    valueRec = bermudanSwaptionRec.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN BDT REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAYER MINUS RECEIVER :", payRec)

    testCases.banner("======= 20% VOLATILITY BERMUDAN SWAPTION EUROPEAN EXERCISE BDT MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.2
    model = FinModelRatesBDT(sigma, numTimeSteps)

    testCases.banner("BDT MODEL BERMUDAN SWAPTION CLASS EUROPEAN EXERCISE")
    valuePay = bermudanSwaptionPay.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN BDT PAY VALUE:", valuePay)

    valueRec = bermudanSwaptionRec.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN BDT REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAYER MINUS RECEIVER :", payRec)
    
    ###########################################################################
    # Now we create the Bermudan swaptions but allow Bermudan exercise
    ###########################################################################

    swapType = FinSwapTypes.PAYER
    exerciseType = FinExerciseTypes.BERMUDAN

    bermudanSwaptionPay = FinLiborBermudanSwaption(settlementDate,
                                                   exerciseDate,
                                                   swapMaturityDate,
                                                   swapType,
                                                   exerciseType,
                                                   swapFixedCoupon,
                                                   swapFixedFrequencyType,
                                                   swapFixedDayCountType)

    swapType = FinSwapTypes.RECEIVER
    bermudanSwaptionRec = FinLiborBermudanSwaption(settlementDate,
                                                   exerciseDate,
                                                   swapMaturityDate,
                                                   swapType,
                                                   exerciseType,
                                                   swapFixedCoupon,
                                                   swapFixedFrequencyType,
                                                   swapFixedDayCountType)

    testCases.banner("======= ZERO VOLATILITY BERMUDAN SWAPTION BERMUDAN EXERCISE BDT MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.000001
    a = 0.01
    model = FinModelRatesBDT(sigma, numTimeSteps)

    testCases.banner("BK MODEL BERMUDAN SWAPTION CLASS BERMUDAN EXERCISE")
    valuePay = bermudanSwaptionPay.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN BDT PAY VALUE:", valuePay)

    valueRec = bermudanSwaptionRec.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN BDT REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAYER MINUS RECEIVER :", payRec)

    testCases.banner("======= 20% VOLATILITY BERMUDAN SWAPTION BERMUDAN EXERCISE BDT MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.20
    a = 0.01
    model = FinModelRatesBDT(sigma, numTimeSteps)

#    print("BDT MODEL BERMUDAN SWAPTION CLASS BERMUDAN EXERCISE")
    valuePay = bermudanSwaptionPay.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN BDT PAY VALUE:", valuePay)

    valueRec = bermudanSwaptionRec.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN BDT REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAYER MINUS RECEIVER :", payRec)
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    # BDT MODEL
    ###########################################################################
    ###########################################################################
    ###########################################################################

    testCases.banner("=======================================================")
    testCases.banner("=======================================================")
    testCases.banner("======================= HW MODEL ======================")
    testCases.banner("=======================================================")
    testCases.banner("=======================================================")

    testCases.banner("====== 0% VOLATILITY EUROPEAN SWAPTION HW MODEL ======")

    sigma = 0.0000001
    a = 0.1
    numTimeSteps = 200
    model = FinModelRatesHW(sigma, a, numTimeSteps)

    valuePay = europeanSwaptionPay.value(valuationDate, liborCurve, model)
    testCases.print("EUROPEAN HW PAY VALUE:", valuePay)

    valueRec = europeanSwaptionRec.value(valuationDate, liborCurve, model)
    testCases.print("EUROPEAN HW REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAYER MINUS RECEIVER :", payRec)

    testCases.banner("===== 20% VOLATILITY EUROPEAN SWAPTION BDT MODEL ======")

    # Used BK with constant short-rate volatility
    sigma = 0.01
    a = 0.01
    model = FinModelRatesHW(sigma, a, numTimeSteps)

    testCases.banner("HW MODEL SWAPTION CLASS EUROPEAN EXERCISE")

    valuePay = europeanSwaptionPay.value(valuationDate, liborCurve, model)
    testCases.print("EUROPEAN HW PAY VALUE:", valuePay)

    valueRec = europeanSwaptionRec.value(valuationDate, liborCurve, model)
    testCases.print("EUROPEAN HW REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAYER MINUS RECEIVER :", payRec)

    ###########################################################################

    # Now we create the Bermudan swaptions but only allow European exercise
    swapType = FinSwapTypes.PAYER
    exerciseType = FinExerciseTypes.EUROPEAN

    bermudanSwaptionPay = FinLiborBermudanSwaption(settlementDate,
                                                   exerciseDate,
                                                   swapMaturityDate,
                                                   swapType,
                                                   exerciseType,
                                                   swapFixedCoupon,
                                                   swapFixedFrequencyType,
                                                   swapFixedDayCountType)

    swapType = FinSwapTypes.RECEIVER
    bermudanSwaptionRec = FinLiborBermudanSwaption(settlementDate,
                                                   exerciseDate,
                                                   swapMaturityDate,
                                                   swapType,
                                                   exerciseType,
                                                   swapFixedCoupon,
                                                   swapFixedFrequencyType,
                                                   swapFixedDayCountType)
   
    testCases.banner("======= 0% VOLATILITY BERMUDAN SWAPTION EUROPEAN EXERCISE HW MODEL ========")

    sigma = 0.000001
    model = FinModelRatesHW(sigma, a, numTimeSteps)

    testCases.banner("BK MODEL BERMUDAN SWAPTION CLASS EUROPEAN EXERCISE")
    valuePay = bermudanSwaptionPay.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN BDT PAY VALUE:", valuePay)

    valueRec = bermudanSwaptionRec.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN BDT REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAYER MINUS RECEIVER :", payRec)

    testCases.banner("======= 100bp VOLATILITY BERMUDAN SWAPTION EUROPEAN EXERCISE HW MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.01
    model = FinModelRatesHW(sigma, a, numTimeSteps)

    testCases.banner("BDT MODEL BERMUDAN SWAPTION CLASS EUROPEAN EXERCISE")
    valuePay = bermudanSwaptionPay.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN BDT PAY VALUE:", valuePay)

    valueRec = bermudanSwaptionRec.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN BDT REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAYER MINUS RECEIVER :", payRec)
    
    ###########################################################################
    # Now we create the Bermudan swaptions but allow Bermudan exercise
    ###########################################################################

    swapType = FinSwapTypes.PAYER
    exerciseType = FinExerciseTypes.BERMUDAN

    bermudanSwaptionPay = FinLiborBermudanSwaption(settlementDate,
                                                   exerciseDate,
                                                   swapMaturityDate,
                                                   swapType,
                                                   exerciseType,
                                                   swapFixedCoupon,
                                                   swapFixedFrequencyType,
                                                   swapFixedDayCountType)

    swapType = FinSwapTypes.RECEIVER
    bermudanSwaptionRec = FinLiborBermudanSwaption(settlementDate,
                                                   exerciseDate,
                                                   swapMaturityDate,
                                                   swapType,
                                                   exerciseType,
                                                   swapFixedCoupon,
                                                   swapFixedFrequencyType,
                                                   swapFixedDayCountType)

    testCases.banner("======= ZERO VOLATILITY BERMUDAN SWAPTION BERMUDAN EXERCISE HW MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.000001
    a = 0.01
    model = FinModelRatesHW(sigma, a, numTimeSteps)

    testCases.banner("HW MODEL BERMUDAN SWAPTION CLASS BERMUDAN EXERCISE")
    valuePay = bermudanSwaptionPay.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN HW PAY VALUE:", valuePay)

    valueRec = bermudanSwaptionRec.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN HW REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAYER MINUS RECEIVER :", payRec)

    testCases.banner("======= 100bps VOLATILITY BERMUDAN SWAPTION BERMUDAN EXERCISE HW MODEL ========")

    # Used BK with constant short-rate volatility
    sigma = 0.01
    a = 0.01
    model = FinModelRatesHW(sigma, a, numTimeSteps)

    testCases.banner("HW MODEL BERMUDAN SWAPTION CLASS BERMUDAN EXERCISE")
    valuePay = bermudanSwaptionPay.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN HW PAY VALUE:", valuePay)

    valueRec = bermudanSwaptionRec.value(valuationDate, liborCurve, model)
    testCases.print("BERMUDAN HW REC VALUE:", valueRec)

    payRec = valuePay - valueRec
    testCases.print("PAYER MINUS RECEIVER :", payRec)
##########################################################################


test_FinLiborBermudanSwaptionBKModel()

testCases.compareTestCases()
