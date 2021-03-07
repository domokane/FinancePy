###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np

import sys
sys.path.append("..")

from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes

from financepy.products.rates.FinIborDeposit import FinIborDeposit
from financepy.products.rates.IborSwap import FinIborSwap
from financepy.products.rates.FinIborSwaption import FinIborSwaption
from financepy.products.rates.FinIborSwaption import FinSwapTypes

from financepy.models.black import FinModelBlack
from financepy.models.black_shifted import FinModelBlackShifted
from financepy.models.sabr import FinModelSABR
from financepy.models.sabr_shifted import FinModelSABRShifted
from financepy.models.rates_hull_white_tree import FinModelRatesHW
from financepy.models.rates_bk_tree import FinModelRatesBK
from financepy.models.rates_bdt_tree import FinModelRatesBDT

from financepy.products.rates.FinIborSingleCurve import IborSingleCurve
from financepy.market.discount.curve_flat import DiscountCurveFlat
from financepy.market.discount.curve_zeros import DiscountCurveZeros
from financepy.market.discount.interpolator import FinInterpTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinIborDepositsAndSwaps(valuation_date):

    depoBasis = DayCountTypes.THIRTY_E_360_ISDA
    depos = []

    spotDays = 0
    settlement_date = valuation_date.addWeekDays(spotDays)
    deposit_rate = 0.05

    depo1 = FinIborDeposit(settlement_date, "1M", deposit_rate, depoBasis)
    depo2 = FinIborDeposit(settlement_date, "3M", deposit_rate, depoBasis)
    depo3 = FinIborDeposit(settlement_date, "6M", deposit_rate, depoBasis)

    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)

    fras = []

    swaps = []
    fixedBasis = DayCountTypes.ACT_365F
    fixedFreq = FrequencyTypes.SEMI_ANNUAL
    fixed_legType = FinSwapTypes.PAY
    
    swap_rate = 0.05
    swap1 = FinIborSwap(settlement_date, "1Y",  fixed_legType, swap_rate, fixedFreq, fixedBasis)
    swap2 = FinIborSwap(settlement_date, "3Y",  fixed_legType, swap_rate, fixedFreq, fixedBasis)
    swap3 = FinIborSwap(settlement_date, "5Y",  fixed_legType, swap_rate, fixedFreq, fixedBasis)

    swaps.append(swap1)
    swaps.append(swap2)
    swaps.append(swap3)

    libor_curve = IborSingleCurve(valuation_date, depos, fras, swaps)

    return libor_curve


##########################################################################

def testFinIborSwaptionModels():

    ##########################################################################
    # COMPARISON OF MODELS
    ##########################################################################

    valuation_date = Date(1, 1, 2011)
    libor_curve = test_FinIborDepositsAndSwaps(valuation_date)

    exerciseDate = Date(1, 1, 2012)
    swapMaturityDate = Date(1, 1, 2017)

    swapFixedFrequencyType = FrequencyTypes.SEMI_ANNUAL
    swapFixedDayCountType = DayCountTypes.ACT_365F

    strikes = np.linspace(0.02, 0.08, 5)

    testCases.header("LAB", "STRIKE", "BLK", "BLK_SHFT", "SABR",
                     "SABR_SHFT", "HW", "BK")

    model1 = FinModelBlack(0.00001)
    model2 = FinModelBlackShifted(0.00001, 0.0)
    model3 = FinModelSABR(0.013, 0.5, 0.5, 0.5)
    model4 = FinModelSABRShifted(0.013, 0.5, 0.5, 0.5, -0.008)
    model5 = FinModelRatesHW(0.00001, 0.00001)
    model6 = FinModelRatesBK(0.01, 0.01)

    settlement_date = valuation_date.addWeekDays(2)

    for k in strikes:
        swaptionType = FinSwapTypes.PAY
        swaption = FinIborSwaption(settlement_date,
                                   exerciseDate,
                                   swapMaturityDate,
                                   swaptionType,
                                   k,
                                   swapFixedFrequencyType,
                                   swapFixedDayCountType)

        swap1 = swaption.value(valuation_date, libor_curve, model1)
        swap2 = swaption.value(valuation_date, libor_curve, model2)
        swap3 = swaption.value(valuation_date, libor_curve, model3)
        swap4 = swaption.value(valuation_date, libor_curve, model4)
        swap5 = swaption.value(valuation_date, libor_curve, model5)
        swap6 = swaption.value(valuation_date, libor_curve, model6)
        testCases.print("PAY", k, swap1, swap2, swap3, swap4, swap5, swap6)

    testCases.header("LABEL", "STRIKE", "BLK", "BLK_SHFTD", "SABR",
                     "SABR_SHFTD", "HW", "BK")

    for k in strikes:
        swaptionType = FinSwapTypes.RECEIVE
        swaption = FinIborSwaption(settlement_date,
                                    exerciseDate,
                                    swapMaturityDate,
                                    swaptionType,
                                    k,
                                    swapFixedFrequencyType,
                                    swapFixedDayCountType)

        swap1 = swaption.value(valuation_date, libor_curve, model1)
        swap2 = swaption.value(valuation_date, libor_curve, model2)
        swap3 = swaption.value(valuation_date, libor_curve, model3)
        swap4 = swaption.value(valuation_date, libor_curve, model4)
        swap5 = swaption.value(valuation_date, libor_curve, model5)
        swap6 = swaption.value(valuation_date, libor_curve, model6)
        testCases.print("REC", k, swap1, swap2, swap3, swap4, swap5, swap6)

###############################################################################


def test_FinIborSwaptionQLExample():

    valuation_date = Date(4, 3, 2014)
    settlement_date = Date(4, 3, 2014)

    depoDCCType = DayCountTypes.THIRTY_E_360_ISDA
    depos = []
    depo = FinIborDeposit(settlement_date, "1W", 0.0023, depoDCCType)
    depos.append(depo)
    depo = FinIborDeposit(settlement_date, "1M", 0.0023, depoDCCType)
    depos.append(depo)
    depo = FinIborDeposit(settlement_date, "3M", 0.0023, depoDCCType)
    depos.append(depo)
    depo = FinIborDeposit(settlement_date, "6M", 0.0023, depoDCCType)
    depos.append(depo)

    # No convexity correction provided so I omit interest rate futures

    swaps = []
    accType = DayCountTypes.ACT_365F
    fixedFreqType = FrequencyTypes.SEMI_ANNUAL
    fixed_legType = FinSwapTypes.PAY
    
    swap = FinIborSwap(settlement_date, "3Y", fixed_legType, 0.00790, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "4Y", fixed_legType, 0.01200, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "5Y", fixed_legType, 0.01570, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "6Y", fixed_legType, 0.01865, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "7Y", fixed_legType, 0.02160, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "8Y", fixed_legType, 0.02350, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "9Y", fixed_legType, 0.02540, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "10Y", fixed_legType, 0.0273, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "15Y", fixed_legType, 0.0297, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "20Y", fixed_legType,  0.0316, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "25Y", fixed_legType, 0.0335, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "30Y", fixed_legType, 0.0354, fixedFreqType, accType)
    swaps.append(swap)

    libor_curve = IborSingleCurve(valuation_date, depos, [], swaps,
                                  FinInterpTypes.LINEAR_ZERO_RATES)

    exerciseDate = settlement_date.addTenor("5Y")
    swapMaturityDate = exerciseDate.addTenor("5Y")
    swapFixedCoupon = 0.040852
    swapFixedFrequencyType = FrequencyTypes.SEMI_ANNUAL
    swapFixedDayCountType = DayCountTypes.THIRTY_E_360_ISDA
    swapFloatFrequencyType = FrequencyTypes.QUARTERLY
    swapFloatDayCountType = DayCountTypes.ACT_360
    swapNotional = 1000000
    swaptionType = FinSwapTypes.PAY

    swaption = FinIborSwaption(settlement_date,
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
    v = swaption.value(settlement_date, libor_curve, model)
    testCases.print(model.__class__, v)

    model = FinModelBlackShifted(0.1533, -0.008)
    v = swaption.value(settlement_date, libor_curve, model)
    testCases.print(model.__class__, v)

    model = FinModelSABR(0.132, 0.5, 0.5, 0.5)
    v = swaption.value(settlement_date, libor_curve, model)
    testCases.print(model.__class__, v)

    model = FinModelSABRShifted(0.352, 0.5, 0.15, 0.15, -0.005)
    v = swaption.value(settlement_date, libor_curve, model)
    testCases.print(model.__class__, v)

    model = FinModelRatesHW(0.010000000, 0.00000000001)
    v = swaption.value(settlement_date, libor_curve, model)
    testCases.print(model.__class__, v)

###############################################################################


def testFinIborCashSettledSwaption():

    testCases.header("LABEL", "VALUE")

    valuation_date = Date(1, 1, 2020)
    settlement_date = Date(1, 1, 2020)

    depoDCCType = DayCountTypes.THIRTY_E_360_ISDA
    depos = []
    depo = FinIborDeposit(settlement_date, "1W", 0.0023, depoDCCType)
    depos.append(depo)
    depo = FinIborDeposit(settlement_date, "1M", 0.0023, depoDCCType)
    depos.append(depo)
    depo = FinIborDeposit(settlement_date, "3M", 0.0023, depoDCCType)
    depos.append(depo)
    depo = FinIborDeposit(settlement_date, "6M", 0.0023, depoDCCType)
    depos.append(depo)

    # No convexity correction provided so I omit interest rate futures

    settlement_date = Date(2, 1, 2020)

    swaps = []
    accType = DayCountTypes.ACT_365F
    fixedFreqType = FrequencyTypes.SEMI_ANNUAL
    fixed_legType = FinSwapTypes.PAY
    
    swap = FinIborSwap(settlement_date, "3Y", fixed_legType, 0.00790, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "4Y", fixed_legType, 0.01200, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "5Y", fixed_legType, 0.01570, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "6Y", fixed_legType, 0.01865, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "7Y", fixed_legType, 0.02160, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "8Y", fixed_legType, 0.02350, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "9Y", fixed_legType, 0.02540, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "10Y", fixed_legType, 0.0273, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "15Y", fixed_legType, 0.0297, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "20Y", fixed_legType, 0.0316, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "25Y", fixed_legType, 0.0335, fixedFreqType, accType)
    swaps.append(swap)
    swap = FinIborSwap(settlement_date, "30Y", fixed_legType, 0.0354, fixedFreqType, accType)
    swaps.append(swap)

    libor_curve = IborSingleCurve(valuation_date, depos, [], swaps,
                                  FinInterpTypes.LINEAR_ZERO_RATES)

    exerciseDate = settlement_date.addTenor("5Y")
    swapMaturityDate = exerciseDate.addTenor("5Y")
    swapFixedCoupon = 0.040852
    swapFixedFrequencyType = FrequencyTypes.SEMI_ANNUAL
    swapFixedDayCountType = DayCountTypes.THIRTY_E_360_ISDA
    swapFloatFrequencyType = FrequencyTypes.QUARTERLY
    swapFloatDayCountType = DayCountTypes.ACT_360
    swapNotional = 1000000
    fixed_legType = FinSwapTypes.PAY

    swaption = FinIborSwaption(settlement_date,
                                exerciseDate,
                                swapMaturityDate,
                                fixed_legType,
                                swapFixedCoupon,
                                swapFixedFrequencyType,
                                swapFixedDayCountType,
                                swapNotional,
                                swapFloatFrequencyType,
                                swapFloatDayCountType)

    model = FinModelBlack(0.1533)
    v = swaption.value(settlement_date, libor_curve, model)
    testCases.print("Swaption No-Arb Value:", v)

    fwdSwapRate1 = libor_curve.swap_rate(exerciseDate,
                                       swapMaturityDate,
                                       swapFixedFrequencyType,
                                       swapFixedDayCountType)

    testCases.print("Curve Fwd Swap Rate:", fwdSwapRate1)

    fwdSwap = FinIborSwap(exerciseDate, 
                          swapMaturityDate, 
                          fixed_legType,
                          swapFixedCoupon, 
                          swapFixedFrequencyType, 
                          swapFixedDayCountType)

    fwdSwapRate2 = fwdSwap.swap_rate(settlement_date, libor_curve)
    testCases.print("Fwd Swap Swap Rate:", fwdSwapRate2)

    model = FinModelBlack(0.1533)

    v = swaption.cashSettledValue(valuation_date,
                                  libor_curve,
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

    valuation_date = Date(1, 1, 2010)
    libor_curve = DiscountCurveFlat(valuation_date, 0.06,
                                    FrequencyTypes.CONTINUOUS,
                                    DayCountTypes.THIRTY_E_360)

    settlement_date = Date(1, 1, 2011)
    exerciseDate = Date(1, 1, 2016)
    maturity_date = Date(1, 1, 2019)

    fixedCoupon = 0.062
    fixedFrequencyType = FrequencyTypes.SEMI_ANNUAL
    fixedDayCountType = DayCountTypes.THIRTY_E_360_ISDA
    notional = 100.0

    # Pricing a PAY
    swaptionType = FinSwapTypes.PAY
    swaption = FinIborSwaption(settlement_date,
                                exerciseDate,
                                maturity_date,
                                swaptionType,
                                fixedCoupon,
                                fixedFrequencyType,
                                fixedDayCountType,
                                notional)

    model = FinModelBlack(0.20)
    v_finpy = swaption.value(valuation_date, libor_curve, model)
    v_matlab = 2.071

    testCases.header("LABEL", "VALUE")
    testCases.print("FP Price:", v_finpy)
    testCases.print("MATLAB Prix:", v_matlab)
    testCases.print("DIFF:", v_finpy - v_matlab)

###############################################################################

    testCases.header("===================================")
    testCases.header("MATLAB EXAMPLE WITH TERM STRUCTURE")
    testCases.header("===================================")

    valuation_date = Date(1, 1, 2010)

    dates = [Date(1, 1, 2011), Date(1, 1, 2012), Date(1, 1, 2013),
             Date(1, 1, 2014), Date(1, 1, 2015)]

    zeroRates = [0.03, 0.034, 0.037, 0.039, 0.040]

    contFreq = FrequencyTypes.CONTINUOUS
    interp_type = FinInterpTypes.LINEAR_ZERO_RATES
    day_count_type = DayCountTypes.THIRTY_E_360

    libor_curve = DiscountCurveZeros(valuation_date, dates,
                                     zeroRates, contFreq,
                                     day_count_type, interp_type)

    settlement_date = Date(1, 1, 2011)
    exerciseDate = Date(1, 1, 2012)
    maturity_date = Date(1, 1, 2017)
    fixedCoupon = 0.03

    fixedFrequencyType = FrequencyTypes.SEMI_ANNUAL
    fixedDayCountType = DayCountTypes.THIRTY_E_360
    floatFrequencyType = FrequencyTypes.SEMI_ANNUAL
    floatDayCountType = DayCountTypes.THIRTY_E_360
    notional = 1000.0

    # Pricing a put
    swaptionType = FinSwapTypes.RECEIVE
    swaption = FinIborSwaption(settlement_date,
                                exerciseDate,
                                maturity_date,
                                swaptionType,
                                fixedCoupon,
                                fixedFrequencyType,
                                fixedDayCountType,
                                notional,
                                floatFrequencyType,
                                floatDayCountType)

    model = FinModelBlack(0.21)
    v_finpy = swaption.value(valuation_date, libor_curve, model)
    v_matlab = 0.5771

    testCases.header("LABEL", "VALUE")
    testCases.print("FP Price:", v_finpy)
    testCases.print("MATLAB Prix:", v_matlab)
    testCases.print("DIFF:", v_finpy - v_matlab)

###############################################################################

    testCases.header("===================================")
    testCases.header("MATLAB EXAMPLE WITH SHIFTED BLACK")
    testCases.header("===================================")

    valuation_date = Date(1, 1, 2016)

    dates = [Date(1, 1, 2017), Date(1, 1, 2018), Date(1, 1, 2019),
             Date(1, 1, 2020), Date(1, 1, 2021)]

    zeroRates = np.array([-0.02, 0.024, 0.047, 0.090, 0.12])/100.0

    contFreq = FrequencyTypes.ANNUAL
    interp_type = FinInterpTypes.LINEAR_ZERO_RATES
    day_count_type = DayCountTypes.THIRTY_E_360

    libor_curve = DiscountCurveZeros(valuation_date, dates, zeroRates,
                                     contFreq, day_count_type, interp_type)

    settlement_date = Date(1, 1, 2016)
    exerciseDate = Date(1, 1, 2017)
    maturity_date = Date(1, 1, 2020)
    fixedCoupon = -0.003

    fixedFrequencyType = FrequencyTypes.SEMI_ANNUAL
    fixedDayCountType = DayCountTypes.THIRTY_E_360_ISDA
    floatFrequencyType = FrequencyTypes.SEMI_ANNUAL
    floatDayCountType = DayCountTypes.THIRTY_E_360_ISDA
    notional = 1000.0

    # Pricing a PAY
    swaptionType = FinSwapTypes.PAY
    swaption = FinIborSwaption(settlement_date,
                                exerciseDate,
                                maturity_date,
                                swaptionType,
                                fixedCoupon,
                                fixedFrequencyType,
                                fixedDayCountType,
                                notional,
                                floatFrequencyType,
                                floatDayCountType)

    model = FinModelBlackShifted(0.31, 0.008)
    v_finpy = swaption.value(valuation_date, libor_curve, model)
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

    valuation_date = Date(1, 1, 2007)

    dates = [Date(1, 1, 2007), Date(1, 7, 2007), Date(1, 1, 2008),
             Date(1, 7, 2008), Date(1, 1, 2009), Date(1, 7, 2009),
             Date(1, 1, 2010), Date(1, 7, 2010),
             Date(1, 1, 2011), Date(1, 7, 2011), Date(1, 1, 2012)]

    zeroRates = np.array([0.075] * 11)
    interp_type = FinInterpTypes.FLAT_FWD_RATES
    day_count_type = DayCountTypes.THIRTY_E_360_ISDA
    contFreq = FrequencyTypes.SEMI_ANNUAL

    libor_curve = DiscountCurveZeros(valuation_date, dates, zeroRates,
                                     contFreq,
                                     day_count_type, interp_type)

    settlement_date = valuation_date
    exerciseDate = Date(1, 1, 2010)
    maturity_date = Date(1, 1, 2012)
    fixedCoupon = 0.04

    fixedFrequencyType = FrequencyTypes.SEMI_ANNUAL
    fixedDayCountType = DayCountTypes.THIRTY_E_360_ISDA
    notional = 100.0

    swaptionType = FinSwapTypes.RECEIVE
    swaption = FinIborSwaption(settlement_date,
                                exerciseDate,
                                maturity_date,
                                swaptionType,
                                fixedCoupon,
                                fixedFrequencyType,
                                fixedDayCountType,
                                notional)

    model = FinModelRatesHW(0.05, 0.01)
    v_finpy = swaption.value(valuation_date, libor_curve, model)
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
    valuation_date = Date(1, 1, 2007)

    dates = [Date(1, 1, 2007), Date(1, 7, 2007), Date(1, 1, 2008),
             Date(1, 7, 2008), Date(1, 1, 2009), Date(1, 7, 2009),
             Date(1, 1, 2010), Date(1, 7, 2010),
             Date(1, 1, 2011), Date(1, 7, 2011), Date(1, 1, 2012)]

    zeroRates = np.array([0.07] * 11)

    interp_type = FinInterpTypes.FLAT_FWD_RATES
    day_count_type = DayCountTypes.THIRTY_E_360_ISDA
    contFreq = FrequencyTypes.SEMI_ANNUAL

    libor_curve = DiscountCurveZeros(valuation_date, dates, zeroRates,
                                     contFreq, day_count_type, interp_type)

    settlement_date = valuation_date
    exerciseDate = Date(1, 1, 2011)
    maturity_date = Date(1, 1, 2012)

    fixedFrequencyType = FrequencyTypes.SEMI_ANNUAL
    fixedDayCountType = DayCountTypes.THIRTY_E_360_ISDA
    notional = 100.0

    model = FinModelRatesBK(0.1, 0.05, 200)

    fixedCoupon = 0.07
    swaptionType = FinSwapTypes.PAY
    swaption = FinIborSwaption(settlement_date,
                                exerciseDate,
                                maturity_date,
                                swaptionType,
                                fixedCoupon,
                                fixedFrequencyType,
                                fixedDayCountType,
                                notional)

    v_finpy = swaption.value(valuation_date, libor_curve, model)
    v_matlab = 0.3634

    testCases.header("LABEL", "VALUE")
    testCases.print("FP Price:", v_finpy)
    testCases.print("MATLAB Prix:", v_matlab)
    testCases.print("DIFF:", v_finpy - v_matlab)

    fixedCoupon = 0.0725
    swaptionType = FinSwapTypes.RECEIVE
    swaption = FinIborSwaption(settlement_date,
                                exerciseDate,
                                maturity_date,
                                swaptionType,
                                fixedCoupon,
                                fixedFrequencyType,
                                fixedDayCountType,
                                notional)

    v_finpy = swaption.value(valuation_date, libor_curve, model)
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

    valuation_date = Date(1, 1, 2007)

    dates = [Date(1, 1, 2007), Date(1, 7, 2007), Date(1, 1, 2008),
             Date(1, 7, 2008), Date(1, 1, 2009), Date(1, 7, 2009),
             Date(1, 1, 2010), Date(1, 7, 2010),
             Date(1, 1, 2011), Date(1, 7, 2011), Date(1, 1, 2012)]

    zeroRates = np.array([0.06] * 11)

    interp_type = FinInterpTypes.FLAT_FWD_RATES
    day_count_type = DayCountTypes.THIRTY_E_360_ISDA
    contFreq = FrequencyTypes.ANNUAL

    libor_curve = DiscountCurveZeros(valuation_date, dates, zeroRates,
                                     contFreq, day_count_type, interp_type)

    settlement_date = valuation_date
    exerciseDate = Date(1, 1, 2012)
    maturity_date = Date(1, 1, 2015)

    fixedFrequencyType = FrequencyTypes.ANNUAL
    fixedDayCountType = DayCountTypes.THIRTY_E_360_ISDA
    notional = 100.0

    fixedCoupon = 0.062
    swaptionType = FinSwapTypes.PAY
    swaption = FinIborSwaption(settlement_date,
                                exerciseDate,
                                maturity_date,
                                swaptionType,
                                fixedCoupon,
                                fixedFrequencyType,
                                fixedDayCountType,
                                notional)

    model = FinModelRatesBDT(0.20, 200)
    v_finpy = swaption.value(valuation_date, libor_curve, model)
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
