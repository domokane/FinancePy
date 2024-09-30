###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys

sys.path.append("..")

import numpy as np
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.rates.ibor_swaption import IborSwaption
from financepy.products.rates.ibor_swaption import SwapTypes
from financepy.models.black import Black
from financepy.models.black_shifted import BlackShifted
from financepy.models.sabr import SABR
from financepy.models.sabr_shifted import SABRShifted
from financepy.models.hw_tree import HWTree
from financepy.models.bk_tree import BKTree
from financepy.models.bdt_tree import BDTTree
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.market.curves.discount_curve_zeros import DiscountCurveZeros
from financepy.market.curves.interpolator import InterpTypes
from FinTestCases import FinTestCases, globalTestCaseMode


test_cases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_ibor_depositsAndSwaps(value_dt):

    depoBasis = DayCountTypes.THIRTY_E_360_ISDA
    depos = []

    spot_days = 0
    settle_dt = value_dt.add_weekdays(spot_days)
    deposit_rate = 0.05

    depo1 = IborDeposit(settle_dt, "1M", deposit_rate, depoBasis)
    depo2 = IborDeposit(settle_dt, "3M", deposit_rate, depoBasis)
    depo3 = IborDeposit(settle_dt, "6M", deposit_rate, depoBasis)

    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)

    fras = []

    swaps = []
    fixedBasis = DayCountTypes.ACT_365F
    fixed_freq = FrequencyTypes.SEMI_ANNUAL
    fixed_leg_type = SwapTypes.PAY

    swap_rate = 0.05
    swap1 = IborSwap(
        settle_dt, "1Y", fixed_leg_type, swap_rate, fixed_freq, fixedBasis
    )
    swap2 = IborSwap(
        settle_dt, "3Y", fixed_leg_type, swap_rate, fixed_freq, fixedBasis
    )
    swap3 = IborSwap(
        settle_dt, "5Y", fixed_leg_type, swap_rate, fixed_freq, fixedBasis
    )

    swaps.append(swap1)
    swaps.append(swap2)
    swaps.append(swap3)

    libor_curve = IborSingleCurve(value_dt, depos, fras, swaps)

    return libor_curve


##########################################################################


def testIborSwaptionModels():

    ##########################################################################
    # COMPARISON OF MODELS
    ##########################################################################

    value_dt = Date(1, 1, 2011)
    libor_curve = test_ibor_depositsAndSwaps(value_dt)

    exercise_dt = Date(1, 1, 2012)
    swap_maturity_dt = Date(1, 1, 2017)

    swap_fixed_freq_type = FrequencyTypes.SEMI_ANNUAL
    swapFixedDayCountType = DayCountTypes.ACT_365F

    strikes = np.linspace(0.02, 0.08, 5)

    test_cases.header(
        "LAB", "STRIKE", "BLK", "BLK_SHFT", "SABR", "SABR_SHFT", "HW", "BK"
    )

    model1 = Black(0.00001)
    model2 = BlackShifted(0.00001, 0.0)
    model3 = SABR(0.013, 0.5, 0.5, 0.5)
    model4 = SABRShifted(0.013, 0.5, 0.5, 0.5, -0.008)
    model5 = HWTree(0.00001, 0.00001)
    model6 = BKTree(0.01, 0.01)

    settle_dt = value_dt.add_weekdays(2)

    for k in strikes:
        swaptionType = SwapTypes.PAY
        swaption = IborSwaption(
            settle_dt,
            exercise_dt,
            swap_maturity_dt,
            swaptionType,
            k,
            swap_fixed_freq_type,
            swapFixedDayCountType,
        )

        swap1 = swaption.value(value_dt, libor_curve, model1)
        swap2 = swaption.value(value_dt, libor_curve, model2)
        swap3 = swaption.value(value_dt, libor_curve, model3)
        swap4 = swaption.value(value_dt, libor_curve, model4)
        swap5 = swaption.value(value_dt, libor_curve, model5)
        swap6 = swaption.value(value_dt, libor_curve, model6)
        test_cases.print("PAY", k, swap1, swap2, swap3, swap4, swap5, swap6)

    test_cases.header(
        "LABEL", "STRIKE", "BLK", "BLK_SHFTD", "SABR", "SABR_SHFTD", "HW", "BK"
    )

    for k in strikes:
        swaptionType = SwapTypes.RECEIVE
        swaption = IborSwaption(
            settle_dt,
            exercise_dt,
            swap_maturity_dt,
            swaptionType,
            k,
            swap_fixed_freq_type,
            swapFixedDayCountType,
        )

        swap1 = swaption.value(value_dt, libor_curve, model1)
        swap2 = swaption.value(value_dt, libor_curve, model2)
        swap3 = swaption.value(value_dt, libor_curve, model3)
        swap4 = swaption.value(value_dt, libor_curve, model4)
        swap5 = swaption.value(value_dt, libor_curve, model5)
        swap6 = swaption.value(value_dt, libor_curve, model6)
        test_cases.print("REC", k, swap1, swap2, swap3, swap4, swap5, swap6)


###############################################################################


def test_IborSwaptionQLExample():

    value_dt = Date(4, 3, 2014)
    settle_dt = Date(4, 3, 2014)

    depoDCCType = DayCountTypes.THIRTY_E_360_ISDA
    depos = []
    depo = IborDeposit(settle_dt, "1W", 0.0023, depoDCCType)
    depos.append(depo)
    depo = IborDeposit(settle_dt, "1M", 0.0023, depoDCCType)
    depos.append(depo)
    depo = IborDeposit(settle_dt, "3M", 0.0023, depoDCCType)
    depos.append(depo)
    depo = IborDeposit(settle_dt, "6M", 0.0023, depoDCCType)
    depos.append(depo)

    # No convexity correction provided so I omit interest rate futures

    swaps = []
    accType = DayCountTypes.ACT_365F
    fixed_freq_type = FrequencyTypes.SEMI_ANNUAL
    fixed_leg_type = SwapTypes.PAY

    swap = IborSwap(
        settle_dt, "3Y", fixed_leg_type, 0.00790, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "4Y", fixed_leg_type, 0.01200, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "5Y", fixed_leg_type, 0.01570, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "6Y", fixed_leg_type, 0.01865, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "7Y", fixed_leg_type, 0.02160, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "8Y", fixed_leg_type, 0.02350, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "9Y", fixed_leg_type, 0.02540, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "10Y", fixed_leg_type, 0.0273, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "15Y", fixed_leg_type, 0.0297, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "20Y", fixed_leg_type, 0.0316, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "25Y", fixed_leg_type, 0.0335, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "30Y", fixed_leg_type, 0.0354, fixed_freq_type, accType
    )
    swaps.append(swap)

    libor_curve = IborSingleCurve(
        value_dt, depos, [], swaps, InterpTypes.LINEAR_ZERO_RATES
    )

    exercise_dt = settle_dt.add_tenor("5Y")
    swap_maturity_dt = exercise_dt.add_tenor("5Y")
    swap_fixed_cpn = 0.040852
    swap_fixed_freq_type = FrequencyTypes.SEMI_ANNUAL
    swapFixedDayCountType = DayCountTypes.THIRTY_E_360_ISDA
    swapFloatFrequencyType = FrequencyTypes.QUARTERLY
    swapFloatDayCountType = DayCountTypes.ACT_360
    swapNotional = 1000000
    swaptionType = SwapTypes.PAY

    swaption = IborSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        swaptionType,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swapFixedDayCountType,
        swapNotional,
        swapFloatFrequencyType,
        swapFloatDayCountType,
    )

    test_cases.header("MODEL", "VALUE")

    model = Black(0.1533)
    v = swaption.value(settle_dt, libor_curve, model)
    test_cases.print(model.__class__, v)

    model = BlackShifted(0.1533, -0.008)
    v = swaption.value(settle_dt, libor_curve, model)
    test_cases.print(model.__class__, v)

    model = SABR(0.132, 0.5, 0.5, 0.5)
    v = swaption.value(settle_dt, libor_curve, model)
    test_cases.print(model.__class__, v)

    model = SABRShifted(0.352, 0.5, 0.15, 0.15, -0.005)
    v = swaption.value(settle_dt, libor_curve, model)
    test_cases.print(model.__class__, v)

    model = HWTree(0.010000000, 0.00000000001)
    v = swaption.value(settle_dt, libor_curve, model)
    test_cases.print(model.__class__, v)


###############################################################################


def testFinIborCashSettledSwaption():

    test_cases.header("LABEL", "VALUE")

    value_dt = Date(1, 1, 2020)
    settle_dt = Date(1, 1, 2020)

    depoDCCType = DayCountTypes.THIRTY_E_360_ISDA
    depos = []
    depo = IborDeposit(settle_dt, "1W", 0.0023, depoDCCType)
    depos.append(depo)
    depo = IborDeposit(settle_dt, "1M", 0.0023, depoDCCType)
    depos.append(depo)
    depo = IborDeposit(settle_dt, "3M", 0.0023, depoDCCType)
    depos.append(depo)
    depo = IborDeposit(settle_dt, "6M", 0.0023, depoDCCType)
    depos.append(depo)

    # No convexity correction provided so I omit interest rate futures

    settle_dt = Date(2, 1, 2020)

    swaps = []
    accType = DayCountTypes.ACT_365F
    fixed_freq_type = FrequencyTypes.SEMI_ANNUAL
    fixed_leg_type = SwapTypes.PAY

    swap = IborSwap(
        settle_dt, "3Y", fixed_leg_type, 0.00790, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "4Y", fixed_leg_type, 0.01200, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "5Y", fixed_leg_type, 0.01570, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "6Y", fixed_leg_type, 0.01865, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "7Y", fixed_leg_type, 0.02160, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "8Y", fixed_leg_type, 0.02350, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "9Y", fixed_leg_type, 0.02540, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "10Y", fixed_leg_type, 0.0273, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "15Y", fixed_leg_type, 0.0297, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "20Y", fixed_leg_type, 0.0316, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "25Y", fixed_leg_type, 0.0335, fixed_freq_type, accType
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt, "30Y", fixed_leg_type, 0.0354, fixed_freq_type, accType
    )
    swaps.append(swap)

    libor_curve = IborSingleCurve(
        value_dt, depos, [], swaps, InterpTypes.LINEAR_ZERO_RATES
    )

    exercise_dt = settle_dt.add_tenor("5Y")
    swap_maturity_dt = exercise_dt.add_tenor("5Y")
    swap_fixed_cpn = 0.040852
    swap_fixed_freq_type = FrequencyTypes.SEMI_ANNUAL
    swapFixedDayCountType = DayCountTypes.THIRTY_E_360_ISDA
    swapFloatFrequencyType = FrequencyTypes.QUARTERLY
    swapFloatDayCountType = DayCountTypes.ACT_360
    swapNotional = 1000000
    fixed_leg_type = SwapTypes.PAY

    swaption = IborSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swapFixedDayCountType,
        swapNotional,
        swapFloatFrequencyType,
        swapFloatDayCountType,
    )

    model = Black(0.1533)
    v = swaption.value(settle_dt, libor_curve, model)
    test_cases.print("Swaption No-Arb Value:", v)

    fwd_swap_rate1 = libor_curve.swap_rate(
        exercise_dt,
        swap_maturity_dt,
        swap_fixed_freq_type,
        swapFixedDayCountType,
    )

    test_cases.print("Curve Fwd Swap Rate:", fwd_swap_rate1)

    fwdSwap = IborSwap(
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swapFixedDayCountType,
    )

    fwd_swap_rate2 = fwdSwap.swap_rate(settle_dt, libor_curve)
    test_cases.print("Fwd Swap Swap Rate:", fwd_swap_rate2)

    model = Black(0.1533)

    v = swaption.cash_settled_value(
        value_dt, libor_curve, fwd_swap_rate2, model
    )

    test_cases.print("Swaption Cash Settled Value:", v)


###############################################################################


def testIborSwaptionMatlabExamples():

    # We value a European swaption using Black's model and try to replicate a
    # ML example at https://fr.mathworks.com/help/fininst/swaptionbyblk.html

    test_cases.header("=======================================")
    test_cases.header("MATLAB EXAMPLE WITH FLAT TERM STRUCTURE")
    test_cases.header("=======================================")

    value_dt = Date(1, 1, 2010)
    libor_curve = DiscountCurveFlat(
        value_dt, 0.06, FrequencyTypes.CONTINUOUS, DayCountTypes.THIRTY_E_360
    )

    settle_dt = Date(1, 1, 2011)
    exercise_dt = Date(1, 1, 2016)
    maturity_dt = Date(1, 1, 2019)

    fixed_cpn = 0.062
    fixed_freq_type = FrequencyTypes.SEMI_ANNUAL
    fixed_dc_type = DayCountTypes.THIRTY_E_360_ISDA
    notional = 100.0

    # Pricing a PAY
    swaptionType = SwapTypes.PAY
    swaption = IborSwaption(
        settle_dt,
        exercise_dt,
        maturity_dt,
        swaptionType,
        fixed_cpn,
        fixed_freq_type,
        fixed_dc_type,
        notional,
    )

    model = Black(0.20)
    v_finpy = swaption.value(value_dt, libor_curve, model)
    v_matlab = 2.071

    test_cases.header("LABEL", "VALUE")
    test_cases.print("FP Price:", v_finpy)
    test_cases.print("MATLAB Prix:", v_matlab)
    test_cases.print("DIFF:", v_finpy - v_matlab)

    ###############################################################################

    test_cases.header("===================================")
    test_cases.header("MATLAB EXAMPLE WITH TERM STRUCTURE")
    test_cases.header("===================================")

    value_dt = Date(1, 1, 2010)

    dates = [
        Date(1, 1, 2011),
        Date(1, 1, 2012),
        Date(1, 1, 2013),
        Date(1, 1, 2014),
        Date(1, 1, 2015),
    ]

    zero_rates = [0.03, 0.034, 0.037, 0.039, 0.040]

    contFreq = FrequencyTypes.CONTINUOUS
    interp_type = InterpTypes.LINEAR_ZERO_RATES
    dc_type = DayCountTypes.THIRTY_E_360

    libor_curve = DiscountCurveZeros(
        value_dt, dates, zero_rates, contFreq, dc_type, interp_type
    )

    settle_dt = Date(1, 1, 2011)
    exercise_dt = Date(1, 1, 2012)
    maturity_dt = Date(1, 1, 2017)
    fixed_cpn = 0.03

    fixed_freq_type = FrequencyTypes.SEMI_ANNUAL
    fixed_dc_type = DayCountTypes.THIRTY_E_360
    float_freq_type = FrequencyTypes.SEMI_ANNUAL
    float_dc_type = DayCountTypes.THIRTY_E_360
    notional = 1000.0

    # Pricing a put
    swaptionType = SwapTypes.RECEIVE
    swaption = IborSwaption(
        settle_dt,
        exercise_dt,
        maturity_dt,
        swaptionType,
        fixed_cpn,
        fixed_freq_type,
        fixed_dc_type,
        notional,
        float_freq_type,
        float_dc_type,
    )

    model = Black(0.21)
    v_finpy = swaption.value(value_dt, libor_curve, model)
    v_matlab = 0.5771

    test_cases.header("LABEL", "VALUE")
    test_cases.print("FP Price:", v_finpy)
    test_cases.print("MATLAB Prix:", v_matlab)
    test_cases.print("DIFF:", v_finpy - v_matlab)

    ###############################################################################

    test_cases.header("===================================")
    test_cases.header("MATLAB EXAMPLE WITH SHIFTED BLACK")
    test_cases.header("===================================")

    value_dt = Date(1, 1, 2016)

    dates = [
        Date(1, 1, 2017),
        Date(1, 1, 2018),
        Date(1, 1, 2019),
        Date(1, 1, 2020),
        Date(1, 1, 2021),
    ]

    zero_rates = np.array([-0.02, 0.024, 0.047, 0.090, 0.12]) / 100.0

    contFreq = FrequencyTypes.ANNUAL
    interp_type = InterpTypes.LINEAR_ZERO_RATES
    dc_type = DayCountTypes.THIRTY_E_360

    libor_curve = DiscountCurveZeros(
        value_dt, dates, zero_rates, contFreq, dc_type, interp_type
    )

    settle_dt = Date(1, 1, 2016)
    exercise_dt = Date(1, 1, 2017)
    maturity_dt = Date(1, 1, 2020)
    fixed_cpn = -0.003

    fixed_freq_type = FrequencyTypes.SEMI_ANNUAL
    fixed_dc_type = DayCountTypes.THIRTY_E_360_ISDA
    float_freq_type = FrequencyTypes.SEMI_ANNUAL
    float_dc_type = DayCountTypes.THIRTY_E_360_ISDA
    notional = 1000.0

    # Pricing a PAY
    swaptionType = SwapTypes.PAY
    swaption = IborSwaption(
        settle_dt,
        exercise_dt,
        maturity_dt,
        swaptionType,
        fixed_cpn,
        fixed_freq_type,
        fixed_dc_type,
        notional,
        float_freq_type,
        float_dc_type,
    )

    model = BlackShifted(0.31, 0.008)
    v_finpy = swaption.value(value_dt, libor_curve, model)
    v_matlab = 12.8301

    test_cases.header("LABEL", "VALUE")
    test_cases.print("FP Price:", v_finpy)
    test_cases.print("MATLAB Prix:", v_matlab)
    test_cases.print("DIFF:", v_finpy - v_matlab)

    ###############################################################################

    test_cases.header("===================================")
    test_cases.header("MATLAB EXAMPLE WITH HULL WHITE")
    test_cases.header("===================================")

    # https://fr.mathworks.com/help/fininst/swaptionbyhw.html

    value_dt = Date(1, 1, 2007)

    dates = [
        Date(1, 1, 2007),
        Date(1, 7, 2007),
        Date(1, 1, 2008),
        Date(1, 7, 2008),
        Date(1, 1, 2009),
        Date(1, 7, 2009),
        Date(1, 1, 2010),
        Date(1, 7, 2010),
        Date(1, 1, 2011),
        Date(1, 7, 2011),
        Date(1, 1, 2012),
    ]

    zero_rates = np.array([0.075] * 11)
    interp_type = InterpTypes.FLAT_FWD_RATES
    dc_type = DayCountTypes.THIRTY_E_360_ISDA
    contFreq = FrequencyTypes.SEMI_ANNUAL

    libor_curve = DiscountCurveZeros(
        value_dt, dates, zero_rates, contFreq, dc_type, interp_type
    )

    settle_dt = value_dt
    exercise_dt = Date(1, 1, 2010)
    maturity_dt = Date(1, 1, 2012)
    fixed_cpn = 0.04

    fixed_freq_type = FrequencyTypes.SEMI_ANNUAL
    fixed_dc_type = DayCountTypes.THIRTY_E_360_ISDA
    notional = 100.0

    swaptionType = SwapTypes.RECEIVE
    swaption = IborSwaption(
        settle_dt,
        exercise_dt,
        maturity_dt,
        swaptionType,
        fixed_cpn,
        fixed_freq_type,
        fixed_dc_type,
        notional,
    )

    model = HWTree(0.05, 0.01)
    v_finpy = swaption.value(value_dt, libor_curve, model)
    v_matlab = 2.9201

    test_cases.header("LABEL", "VALUE")
    test_cases.print("FP Price:", v_finpy)
    test_cases.print("MATLAB Prix:", v_matlab)
    test_cases.print("DIFF:", v_finpy - v_matlab)

    ###############################################################################

    test_cases.header("====================================")
    test_cases.header("MATLAB EXAMPLE WITH BLACK KARASINSKI")
    test_cases.header("====================================")

    # https://fr.mathworks.com/help/fininst/swaptionbybk.html
    value_dt = Date(1, 1, 2007)

    dates = [
        Date(1, 1, 2007),
        Date(1, 7, 2007),
        Date(1, 1, 2008),
        Date(1, 7, 2008),
        Date(1, 1, 2009),
        Date(1, 7, 2009),
        Date(1, 1, 2010),
        Date(1, 7, 2010),
        Date(1, 1, 2011),
        Date(1, 7, 2011),
        Date(1, 1, 2012),
    ]

    zero_rates = np.array([0.07] * 11)

    interp_type = InterpTypes.FLAT_FWD_RATES
    dc_type = DayCountTypes.THIRTY_E_360_ISDA
    contFreq = FrequencyTypes.SEMI_ANNUAL

    libor_curve = DiscountCurveZeros(
        value_dt, dates, zero_rates, contFreq, dc_type, interp_type
    )

    settle_dt = value_dt
    exercise_dt = Date(1, 1, 2011)
    maturity_dt = Date(1, 1, 2012)

    fixed_freq_type = FrequencyTypes.SEMI_ANNUAL
    fixed_dc_type = DayCountTypes.THIRTY_E_360_ISDA
    notional = 100.0

    model = BKTree(0.1, 0.05, 200)

    fixed_cpn = 0.07
    swaptionType = SwapTypes.PAY
    swaption = IborSwaption(
        settle_dt,
        exercise_dt,
        maturity_dt,
        swaptionType,
        fixed_cpn,
        fixed_freq_type,
        fixed_dc_type,
        notional,
    )

    v_finpy = swaption.value(value_dt, libor_curve, model)
    v_matlab = 0.3634

    test_cases.header("LABEL", "VALUE")
    test_cases.print("FP Price:", v_finpy)
    test_cases.print("MATLAB Prix:", v_matlab)
    test_cases.print("DIFF:", v_finpy - v_matlab)

    fixed_cpn = 0.0725
    swaptionType = SwapTypes.RECEIVE
    swaption = IborSwaption(
        settle_dt,
        exercise_dt,
        maturity_dt,
        swaptionType,
        fixed_cpn,
        fixed_freq_type,
        fixed_dc_type,
        notional,
    )

    v_finpy = swaption.value(value_dt, libor_curve, model)
    v_matlab = 0.4798

    test_cases.header("LABEL", "VALUE")
    test_cases.print("FP Price:", v_finpy)
    test_cases.print("MATLAB Prix:", v_matlab)
    test_cases.print("DIFF:", v_finpy - v_matlab)

    ###############################################################################

    test_cases.header("====================================")
    test_cases.header("MATLAB EXAMPLE WITH BLACK-DERMAN-TOY")
    test_cases.header("====================================")

    # https://fr.mathworks.com/help/fininst/swaptionbybdt.html

    value_dt = Date(1, 1, 2007)

    dates = [
        Date(1, 1, 2007),
        Date(1, 7, 2007),
        Date(1, 1, 2008),
        Date(1, 7, 2008),
        Date(1, 1, 2009),
        Date(1, 7, 2009),
        Date(1, 1, 2010),
        Date(1, 7, 2010),
        Date(1, 1, 2011),
        Date(1, 7, 2011),
        Date(1, 1, 2012),
    ]

    zero_rates = np.array([0.06] * 11)

    interp_type = InterpTypes.FLAT_FWD_RATES
    dc_type = DayCountTypes.THIRTY_E_360_ISDA
    contFreq = FrequencyTypes.ANNUAL

    libor_curve = DiscountCurveZeros(
        value_dt, dates, zero_rates, contFreq, dc_type, interp_type
    )

    settle_dt = value_dt
    exercise_dt = Date(1, 1, 2012)
    maturity_dt = Date(1, 1, 2015)

    fixed_freq_type = FrequencyTypes.ANNUAL
    fixed_dc_type = DayCountTypes.THIRTY_E_360_ISDA
    notional = 100.0

    fixed_cpn = 0.062
    swaptionType = SwapTypes.PAY
    swaption = IborSwaption(
        settle_dt,
        exercise_dt,
        maturity_dt,
        swaptionType,
        fixed_cpn,
        fixed_freq_type,
        fixed_dc_type,
        notional,
    )

    model = BDTTree(0.20, 200)
    v_finpy = swaption.value(value_dt, libor_curve, model)
    v_matlab = 2.0592

    test_cases.header("LABEL", "VALUE")
    test_cases.print("FP Price:", v_finpy)
    test_cases.print("MATLAB Prix:", v_matlab)
    test_cases.print("DIFF:", v_finpy - v_matlab)


###############################################################################


testIborSwaptionModels()
testFinIborCashSettledSwaption()
testIborSwaptionMatlabExamples()
test_IborSwaptionQLExample()

test_cases.compareTestCases()
