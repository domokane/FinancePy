###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np

import sys
sys.path.append("..")

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes

from financepy.products.rates.FinIborDeposit import FinIborDeposit
from financepy.products.rates.FinIborSwap import FinIborSwap
from financepy.products.rates.FinIborSwaption import FinIborSwaption
from financepy.products.rates.FinIborSwaption import FinSwapTypes

from financepy.models.FinModelBlack import FinModelBlack
from financepy.models.FinModelBlackShifted import FinModelBlackShifted
from financepy.models.FinModelSABR import FinModelSABR
from financepy.models.FinModelSABRShifted import FinModelSABRShifted
from financepy.models.FinModelRatesHW import FinModelRatesHW
from financepy.models.FinModelRatesBK import FinModelRatesBK
from financepy.models.FinModelRatesBDT import FinModelRatesBDT

from financepy.products.rates.FinIborSingleCurve import FinIborSingleCurve
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.market.curves.FinDiscountCurveZeros import FinDiscountCurveZeros
from financepy.market.curves.FinInterpolator import FinInterpTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinIborDepositsAndSwaps(valuationDate):

    depoBasis = FinDayCountTypes.THIRTY_E_360_ISDA
    depos = []

    spotDays = 0
    settlementDate = valuationDate.addWeekDays(spotDays)
    depositRate = 0.05

    depo1 = FinIborDeposit(settlementDate, "1M", depositRate, depoBasis)
    depo2 = FinIborDeposit(settlementDate, "3M", depositRate, depoBasis)
    depo3 = FinIborDeposit(settlementDate, "6M", depositRate, depoBasis)

    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)

    fras = []

    swaps = []
    fixedBasis = FinDayCountTypes.ACT_365F
    fixedFreq = FinFrequencyTypes.SEMI_ANNUAL
    fixedLegType = FinSwapTypes.PAY
    
    swapRate = 0.05
    swap1 = FinIborSwap(settlementDate, "1Y",  fixedLegType, swapRate, fixedFreq, fixedBasis)
    swap2 = FinIborSwap(settlementDate, "3Y",  fixedLegType, swapRate, fixedFreq, fixedBasis)
    swap3 = FinIborSwap(settlementDate, "5Y",  fixedLegType, swapRate, fixedFreq, fixedBasis)

    swaps.append(swap1)
    swaps.append(swap2)
    swaps.append(swap3)

    liborCurve = FinIborSingleCurve(valuationDate, depos, fras, swaps)

    return liborCurve


##########################################################################

def testFinIborSwaptionModels():

    ##########################################################################
    # COMPARISON OF MODELS
    ##########################################################################

    valuationDate = FinDate(1, 1, 2011)
    liborCurve = test_FinIborDepositsAndSwaps(valuationDate)

    exerciseDate = FinDate(1, 1, 2012)
    swapMaturityDate = FinDate(1, 1, 2017)

    swapFixedFrequencyType = FinFrequencyTypes.SEMI_ANNUAL
    swapFixedDayCountType = FinDayCountTypes.ACT_365F

    strikes = np.linspace(0.02, 0.08, 5)

    testCases.header("LAB", "STRIKE", "BLK", "BLK_SHFT", "SABR",
                     "SABR_SHFT", "HW", "BK")

    model1 = FinModelBlack(0.00001)
    model2 = FinModelBlackShifted(0.00001, 0.0)
    model3 = FinModelSABR(0.013, 0.5, 0.5, 0.5)
    model4 = FinModelSABRShifted(0.013, 0.5, 0.5, 0.5, -0.008)
    model5 = FinModelRatesHW(0.00001, 0.00001)
    model6 = FinModelRatesBK(0.01, 0.01)

    settlementDate = valuationDate.addWeekDays(2)

    for k in strikes:
        swaptionType = FinSwapTypes.PAY
        swaption = FinIborSwaption(settlementDate,
                                   exerciseDate,
                                   swapMaturityDate,
                                   swaptionType,
                                   k,
                                   swapFixedFrequencyType,
                                   swapFixedDayCountType)

        swap1 = swaption.value(valuationDate, liborCurve, model1)
        swap2 = swaption.value(valuationDate, liborCurve, model2)
        swap3 = swaption.value(valuationDate, liborCurve, model3)
        swap4 = swaption.value(valuationDate, liborCurve, model4)
        swap5 = swaption.value(valuationDate, liborCurve, model5)
        swap6 = swaption.value(valuationDate, liborCurve, model6)
        testCases.print("PAY", k, swap1, swap2, swap3, swap4, swap5, swap6)

    testCases.header("LABEL", "STRIKE", "BLK", "BLK_SHFTD", "SABR",
                     "SABR_SHFTD", "HW", "BK")

    for k in strikes:
        swaptionType = FinSwapTypes.RECEIVE
        swaption = FinIborSwaption(settlementDate,
                                    exerciseDate,
                                    swapMaturityDate,
                                    swaptionType,
                                    k,
                                    swapFixedFrequencyType,
                                    swapFixedDayCountType)

        swap1 = swaption.value(valuationDate, liborCurve, model1)
        swap2 = swaption.value(valuationDate, liborCurve, model2)
        swap3 = swaption.value(valuationDate, liborCurve, model3)
        swap4 = swaption.value(valuationDate, liborCurve, model4)
        swap5 = swaption.value(valuationDate, liborCurve, model5)
        swap6 = swaption.value(valuationDate, liborCurve, model6)
        testCases.print("REC", k, swap1, swap2, swap3, swap4, swap5, swap6)

###############################################################################


def test_FinIborSwaptionQLExample():

    valuationDate = FinDate(4, 3, 2014)
    settlementDate = FinDate(4, 3, 2014)

    depoDCCType = FinDayCountTypes.THIRTY_E_360_ISDA
    depos = []
    depo = FinIborDeposit(settlementDate, "1W", 0.0023, depoDCCType)
    depos.append(depo)
    depo = FinIborDeposit(settlementDate, "1M", 0.0023, depoDCCType)
    depos.append(depo)
    depo = FinIborDeposit(settlementDate, "3M", 0.0023, depoDCCType)
    depos.append(depo)
    depo = FinIborDeposit(settlementDate, "6M", 0.0023, depoDCCType)
    depos.append(depo)

    # No convexity correction provided so I omit interest rate futures

    swaps = []
    accType = FinDayCountTypes.ACT_365F
    fixedFreqType = FinFrequencyTypes.SEMI_ANNUAL
    fixedLegType = FinSwapTypes.PAY
    
    swap = FinIborSwap(settlementDate, "3Y", fixedLegType, 0.00790, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "4Y", fixedLegType, 0.01200, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "5Y", fixedLegType, 0.01570, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "6Y", fixedLegType, 0.01865, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "7Y", fixedLegType, 0.02160, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "8Y", fixedLegType, 0.02350, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "9Y", fixedLegType, 0.02540, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "10Y", fixedLegType, 0.0273, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "15Y", fixedLegType, 0.0297, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "20Y", fixedLegType,  0.0316, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "25Y", fixedLegType, 0.0335, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "30Y", fixedLegType, 0.0354, fixedFreqType, accType)
    swaps.append(swap)

    liborCurve = FinIborSingleCurve(valuationDate, depos, [], swaps,
                                    FinInterpTypes.LINEAR_ZERO_RATES)

    exerciseDate = settlementDate.addTenor("5Y")
    swapMaturityDate = exerciseDate.addTenor("5Y")
    swapFixedCoupon = 0.040852
    swapFixedFrequencyType = FinFrequencyTypes.SEMI_ANNUAL
    swapFixedDayCountType = FinDayCountTypes.THIRTY_E_360_ISDA
    swapFloatFrequencyType = FinFrequencyTypes.QUARTERLY
    swapFloatDayCountType = FinDayCountTypes.ACT_360
    swapNotional = 1000000
    swaptionType = FinSwapTypes.PAY

    swaption = FinIborSwaption(settlementDate,
                               exerciseDate,
                               swapMaturityDate,
                               swaptionType,
                               swapFixedCoupon,
                               swapFixedFrequencyType,
                               swapFixedDayCountType,
                               swapNotional,
                               swapFloatFrequencyType,
                               swapFloatDayCountType)

    testCases.header("MODEL", "VALUE")

    model = FinModelBlack(0.1533)
    v = swaption.value(settlementDate, liborCurve, model)
    testCases.print(model.__class__, v)

    model = FinModelBlackShifted(0.1533, -0.008)
    v = swaption.value(settlementDate, liborCurve, model)
    testCases.print(model.__class__, v)

    model = FinModelSABR(0.132, 0.5, 0.5, 0.5)
    v = swaption.value(settlementDate, liborCurve, model)
    testCases.print(model.__class__, v)

    model = FinModelSABRShifted(0.352, 0.5, 0.15, 0.15, -0.005)
    v = swaption.value(settlementDate, liborCurve, model)
    testCases.print(model.__class__, v)

    model = FinModelRatesHW(0.010000000, 0.00000000001)
    v = swaption.value(settlementDate, liborCurve, model)
    testCases.print(model.__class__, v)

###############################################################################


def testFinIborCashSettledSwaption():

    testCases.header("LABEL", "VALUE")

    valuationDate = FinDate(1, 1, 2020)
    settlementDate = FinDate(1, 1, 2020)

    depoDCCType = FinDayCountTypes.THIRTY_E_360_ISDA
    depos = []
    depo = FinIborDeposit(settlementDate, "1W", 0.0023, depoDCCType)
    depos.append(depo)
    depo = FinIborDeposit(settlementDate, "1M", 0.0023, depoDCCType)
    depos.append(depo)
    depo = FinIborDeposit(settlementDate, "3M", 0.0023, depoDCCType)
    depos.append(depo)
    depo = FinIborDeposit(settlementDate, "6M", 0.0023, depoDCCType)
    depos.append(depo)

    # No convexity correction provided so I omit interest rate futures

    settlementDate = FinDate(2, 1, 2020)

    swaps = []
    accType = FinDayCountTypes.ACT_365F
    fixedFreqType = FinFrequencyTypes.SEMI_ANNUAL
    fixedLegType = FinSwapTypes.PAY
    
    swap = FinIborSwap(settlementDate, "3Y", fixedLegType, 0.00790, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "4Y", fixedLegType, 0.01200, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "5Y", fixedLegType, 0.01570, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "6Y", fixedLegType, 0.01865, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "7Y", fixedLegType, 0.02160, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "8Y", fixedLegType, 0.02350, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "9Y", fixedLegType, 0.02540, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "10Y", fixedLegType, 0.0273, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "15Y", fixedLegType, 0.0297, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "20Y", fixedLegType, 0.0316, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "25Y", fixedLegType, 0.0335, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlementDate, "30Y", fixedLegType, 0.0354, fixedFreqType, accType)
    swaps.append(swap)

    liborCurve = FinIborSingleCurve(valuationDate, depos, [], swaps,
                               FinInterpTypes.LINEAR_ZERO_RATES)

    exerciseDate = settlementDate.addTenor("5Y")
    swapMaturityDate = exerciseDate.addTenor("5Y")
    swapFixedCoupon = 0.040852
    swapFixedFrequencyType = FinFrequencyTypes.SEMI_ANNUAL
    swapFixedDayCountType = FinDayCountTypes.THIRTY_E_360_ISDA
    swapFloatFrequencyType = FinFrequencyTypes.QUARTERLY
    swapFloatDayCountType = FinDayCountTypes.ACT_360
    swapNotional = 1000000
    fixedLegType = FinSwapTypes.PAY

    swaption = FinIborSwaption(settlementDate,
                                exerciseDate,
                                swapMaturityDate,
                                fixedLegType,
                                swapFixedCoupon,
                                swapFixedFrequencyType,
                                swapFixedDayCountType,
                                swapNotional,
                                swapFloatFrequencyType,
                                swapFloatDayCountType)

    model = FinModelBlack(0.1533)
    v = swaption.value(settlementDate, liborCurve, model)
    testCases.print("Swaption No-Arb Value:", v)

    fwdSwapRate1 = liborCurve.swapRate(exerciseDate,
                                       swapMaturityDate,
                                       swapFixedFrequencyType,
                                       swapFixedDayCountType)

    testCases.print("Curve Fwd Swap Rate:", fwdSwapRate1)

    fwdSwap = FinIborSwap(exerciseDate, 
                          swapMaturityDate, 
                          fixedLegType, 
                          swapFixedCoupon, 
                          swapFixedFrequencyType, 
                          swapFixedDayCountType)

    fwdSwapRate2 = fwdSwap.swapRate(settlementDate, liborCurve)
    testCases.print("Fwd Swap Swap Rate:", fwdSwapRate2)

    model = FinModelBlack(0.1533)

    v = swaption.cashSettledValue(valuationDate,
                                  liborCurve,
                                  fwdSwapRate2,
                                  model)

    testCases.print("Swaption Cash Settled Value:", v)

###############################################################################


def testFinIborSwaptionMatlabExamples():

    # We value a European swaption using Black's model and try to replicate a
    # ML example at https://fr.mathworks.com/help/fininst/swaptionbyblk.html

    testCases.header("=======================================")
    testCases.header("MATLAB EXAMPLE WITH FLAT TERM STRUCTURE")
    testCases.header("=======================================")

    valuationDate = FinDate(1, 1, 2010)
    liborCurve = FinDiscountCurveFlat(valuationDate, 0.06,
                                      FinFrequencyTypes.CONTINUOUS,
                                      FinDayCountTypes.THIRTY_E_360)

    settlementDate = FinDate(1, 1, 2011)
    exerciseDate = FinDate(1, 1, 2016)
    maturityDate = FinDate(1, 1, 2019)

    fixedCoupon = 0.062
    fixedFrequencyType = FinFrequencyTypes.SEMI_ANNUAL
    fixedDayCountType = FinDayCountTypes.THIRTY_E_360_ISDA
    notional = 100.0

    # Pricing a PAY
    swaptionType = FinSwapTypes.PAY
    swaption = FinIborSwaption(settlementDate,
                                exerciseDate,
                                maturityDate,
                                swaptionType,
                                fixedCoupon,
                                fixedFrequencyType,
                                fixedDayCountType,
                                notional)

    model = FinModelBlack(0.20)
    v_finpy = swaption.value(valuationDate, liborCurve, model)
    v_matlab = 2.071

    testCases.header("LABEL", "VALUE")
    testCases.print("FP Price:", v_finpy)
    testCases.print("MATLAB Prix:", v_matlab)
    testCases.print("DIFF:", v_finpy - v_matlab)

###############################################################################

    testCases.header("===================================")
    testCases.header("MATLAB EXAMPLE WITH TERM STRUCTURE")
    testCases.header("===================================")

    valuationDate = FinDate(1, 1, 2010)

    dates = [FinDate(1, 1, 2011), FinDate(1, 1, 2012), FinDate(1, 1, 2013),
             FinDate(1, 1, 2014), FinDate(1, 1, 2015)]

    zeroRates = [0.03, 0.034, 0.037, 0.039, 0.040]

    contFreq = FinFrequencyTypes.CONTINUOUS
    interpType = FinInterpTypes.LINEAR_ZERO_RATES
    dayCountType = FinDayCountTypes.THIRTY_E_360

    liborCurve = FinDiscountCurveZeros(valuationDate, dates,
                                       zeroRates, contFreq,
                                       dayCountType, interpType)

    settlementDate = FinDate(1, 1, 2011)
    exerciseDate = FinDate(1, 1, 2012)
    maturityDate = FinDate(1, 1, 2017)
    fixedCoupon = 0.03

    fixedFrequencyType = FinFrequencyTypes.SEMI_ANNUAL
    fixedDayCountType = FinDayCountTypes.THIRTY_E_360
    floatFrequencyType = FinFrequencyTypes.SEMI_ANNUAL
    floatDayCountType = FinDayCountTypes.THIRTY_E_360
    notional = 1000.0

    # Pricing a put
    swaptionType = FinSwapTypes.RECEIVE
    swaption = FinIborSwaption(settlementDate,
                                exerciseDate,
                                maturityDate,
                                swaptionType,
                                fixedCoupon,
                                fixedFrequencyType,
                                fixedDayCountType,
                                notional,
                                floatFrequencyType,
                                floatDayCountType)

    model = FinModelBlack(0.21)
    v_finpy = swaption.value(valuationDate, liborCurve, model)
    v_matlab = 0.5771

    testCases.header("LABEL", "VALUE")
    testCases.print("FP Price:", v_finpy)
    testCases.print("MATLAB Prix:", v_matlab)
    testCases.print("DIFF:", v_finpy - v_matlab)

###############################################################################

    testCases.header("===================================")
    testCases.header("MATLAB EXAMPLE WITH SHIFTED BLACK")
    testCases.header("===================================")

    valuationDate = FinDate(1, 1, 2016)

    dates = [FinDate(1, 1, 2017), FinDate(1, 1, 2018), FinDate(1, 1, 2019),
             FinDate(1, 1, 2020), FinDate(1, 1, 2021)]

    zeroRates = np.array([-0.02, 0.024, 0.047, 0.090, 0.12])/100.0

    contFreq = FinFrequencyTypes.ANNUAL
    interpType = FinInterpTypes.LINEAR_ZERO_RATES
    dayCountType = FinDayCountTypes.THIRTY_E_360

    liborCurve = FinDiscountCurveZeros(valuationDate, dates, zeroRates,
                                       contFreq, dayCountType, interpType)

    settlementDate = FinDate(1, 1, 2016)
    exerciseDate = FinDate(1, 1, 2017)
    maturityDate = FinDate(1, 1, 2020)
    fixedCoupon = -0.003

    fixedFrequencyType = FinFrequencyTypes.SEMI_ANNUAL
    fixedDayCountType = FinDayCountTypes.THIRTY_E_360_ISDA
    floatFrequencyType = FinFrequencyTypes.SEMI_ANNUAL
    floatDayCountType = FinDayCountTypes.THIRTY_E_360_ISDA
    notional = 1000.0

    # Pricing a PAY
    swaptionType = FinSwapTypes.PAY
    swaption = FinIborSwaption(settlementDate,
                                exerciseDate,
                                maturityDate,
                                swaptionType,
                                fixedCoupon,
                                fixedFrequencyType,
                                fixedDayCountType,
                                notional,
                                floatFrequencyType,
                                floatDayCountType)

    model = FinModelBlackShifted(0.31, 0.008)
    v_finpy = swaption.value(valuationDate, liborCurve, model)
    v_matlab = 12.8301

    testCases.header("LABEL", "VALUE")
    testCases.print("FP Price:", v_finpy)
    testCases.print("MATLAB Prix:", v_matlab)
    testCases.print("DIFF:", v_finpy - v_matlab)

###############################################################################

    testCases.header("===================================")
    testCases.header("MATLAB EXAMPLE WITH HULL WHITE")
    testCases.header("===================================")

    # https://fr.mathworks.com/help/fininst/swaptionbyhw.html

    valuationDate = FinDate(1, 1, 2007)

    dates = [FinDate(1, 1, 2007), FinDate(1, 7, 2007), FinDate(1, 1, 2008),
             FinDate(1, 7, 2008), FinDate(1, 1, 2009), FinDate(1, 7, 2009),
             FinDate(1, 1, 2010), FinDate(1, 7, 2010),
             FinDate(1, 1, 2011), FinDate(1, 7, 2011), FinDate(1, 1, 2012)]

    zeroRates = np.array([0.075] * 11)
    interpType = FinInterpTypes.FLAT_FWD_RATES
    dayCountType = FinDayCountTypes.THIRTY_E_360_ISDA
    contFreq = FinFrequencyTypes.SEMI_ANNUAL

    liborCurve = FinDiscountCurveZeros(valuationDate, dates, zeroRates,
                                       contFreq,
                                       dayCountType, interpType)

    settlementDate = valuationDate
    exerciseDate = FinDate(1, 1, 2010)
    maturityDate = FinDate(1, 1, 2012)
    fixedCoupon = 0.04

    fixedFrequencyType = FinFrequencyTypes.SEMI_ANNUAL
    fixedDayCountType = FinDayCountTypes.THIRTY_E_360_ISDA
    notional = 100.0

    swaptionType = FinSwapTypes.RECEIVE
    swaption = FinIborSwaption(settlementDate,
                                exerciseDate,
                                maturityDate,
                                swaptionType,
                                fixedCoupon,
                                fixedFrequencyType,
                                fixedDayCountType,
                                notional)

    model = FinModelRatesHW(0.05, 0.01)
    v_finpy = swaption.value(valuationDate, liborCurve, model)
    v_matlab = 2.9201

    testCases.header("LABEL", "VALUE")
    testCases.print("FP Price:", v_finpy)
    testCases.print("MATLAB Prix:", v_matlab)
    testCases.print("DIFF:", v_finpy - v_matlab)

###############################################################################

    testCases.header("====================================")
    testCases.header("MATLAB EXAMPLE WITH BLACK KARASINSKI")
    testCases.header("====================================")

    # https://fr.mathworks.com/help/fininst/swaptionbybk.html
    valuationDate = FinDate(1, 1, 2007)

    dates = [FinDate(1, 1, 2007), FinDate(1, 7, 2007), FinDate(1, 1, 2008),
             FinDate(1, 7, 2008), FinDate(1, 1, 2009), FinDate(1, 7, 2009),
             FinDate(1, 1, 2010), FinDate(1, 7, 2010),
             FinDate(1, 1, 2011), FinDate(1, 7, 2011), FinDate(1, 1, 2012)]

    zeroRates = np.array([0.07] * 11)

    interpType = FinInterpTypes.FLAT_FWD_RATES
    dayCountType = FinDayCountTypes.THIRTY_E_360_ISDA
    contFreq = FinFrequencyTypes.SEMI_ANNUAL

    liborCurve = FinDiscountCurveZeros(valuationDate, dates, zeroRates,
                                       contFreq, dayCountType, interpType)

    settlementDate = valuationDate
    exerciseDate = FinDate(1, 1, 2011)
    maturityDate = FinDate(1, 1, 2012)

    fixedFrequencyType = FinFrequencyTypes.SEMI_ANNUAL
    fixedDayCountType = FinDayCountTypes.THIRTY_E_360_ISDA
    notional = 100.0

    model = FinModelRatesBK(0.1, 0.05, 200)

    fixedCoupon = 0.07
    swaptionType = FinSwapTypes.PAY
    swaption = FinIborSwaption(settlementDate,
                                exerciseDate,
                                maturityDate,
                                swaptionType,
                                fixedCoupon,
                                fixedFrequencyType,
                                fixedDayCountType,
                                notional)

    v_finpy = swaption.value(valuationDate, liborCurve, model)
    v_matlab = 0.3634

    testCases.header("LABEL", "VALUE")
    testCases.print("FP Price:", v_finpy)
    testCases.print("MATLAB Prix:", v_matlab)
    testCases.print("DIFF:", v_finpy - v_matlab)

    fixedCoupon = 0.0725
    swaptionType = FinSwapTypes.RECEIVE
    swaption = FinIborSwaption(settlementDate,
                                exerciseDate,
                                maturityDate,
                                swaptionType,
                                fixedCoupon,
                                fixedFrequencyType,
                                fixedDayCountType,
                                notional)

    v_finpy = swaption.value(valuationDate, liborCurve, model)
    v_matlab = 0.4798

    testCases.header("LABEL", "VALUE")
    testCases.print("FP Price:", v_finpy)
    testCases.print("MATLAB Prix:", v_matlab)
    testCases.print("DIFF:", v_finpy - v_matlab)

###############################################################################

    testCases.header("====================================")
    testCases.header("MATLAB EXAMPLE WITH BLACK-DERMAN-TOY")
    testCases.header("====================================")

    # https://fr.mathworks.com/help/fininst/swaptionbybdt.html

    valuationDate = FinDate(1, 1, 2007)

    dates = [FinDate(1, 1, 2007), FinDate(1, 7, 2007), FinDate(1, 1, 2008),
             FinDate(1, 7, 2008), FinDate(1, 1, 2009), FinDate(1, 7, 2009),
             FinDate(1, 1, 2010), FinDate(1, 7, 2010),
             FinDate(1, 1, 2011), FinDate(1, 7, 2011), FinDate(1, 1, 2012)]

    zeroRates = np.array([0.06] * 11)

    interpType = FinInterpTypes.FLAT_FWD_RATES
    dayCountType = FinDayCountTypes.THIRTY_E_360_ISDA
    contFreq = FinFrequencyTypes.ANNUAL

    liborCurve = FinDiscountCurveZeros(valuationDate, dates, zeroRates,
                                       contFreq, dayCountType, interpType)

    settlementDate = valuationDate
    exerciseDate = FinDate(1, 1, 2012)
    maturityDate = FinDate(1, 1, 2015)

    fixedFrequencyType = FinFrequencyTypes.ANNUAL
    fixedDayCountType = FinDayCountTypes.THIRTY_E_360_ISDA
    notional = 100.0

    fixedCoupon = 0.062
    swaptionType = FinSwapTypes.PAY
    swaption = FinIborSwaption(settlementDate,
                                exerciseDate,
                                maturityDate,
                                swaptionType,
                                fixedCoupon,
                                fixedFrequencyType,
                                fixedDayCountType,
                                notional)

    model = FinModelRatesBDT(0.20, 200)
    v_finpy = swaption.value(valuationDate, liborCurve, model)
    v_matlab = 2.0592

    testCases.header("LABEL", "VALUE")
    testCases.print("FP Price:", v_finpy)
    testCases.print("MATLAB Prix:", v_matlab)
    testCases.print("DIFF:", v_finpy - v_matlab)


###############################################################################


testFinIborSwaptionModels()
testFinIborCashSettledSwaption()
testFinIborSwaptionMatlabExamples()
test_FinIborSwaptionQLExample()

testCases.compareTestCases()
