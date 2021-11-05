###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.bdt_tree import BDTTree
from financepy.models.hw_tree import HWTree
from financepy.models.bk_tree import BKTree
from financepy.models.black import Black
from financepy.products.rates.bermudan_swaption import IborBermudanSwaption
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.rates.ibor_swaption import IborSwaption
from financepy.utils.global_types import FinExerciseTypes
from financepy.utils.global_types import SwapTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date


valuation_date = Date(1, 1, 2011)
settlement_date = valuation_date
exercise_date = settlement_date.add_years(1)
swapMaturityDate = settlement_date.add_years(4)

swapFixedCoupon = 0.060
swapFixedFrequencyType = FrequencyTypes.SEMI_ANNUAL
swapFixedDayCountType = DayCountTypes.ACT_365F

libor_curve = DiscountCurveFlat(valuation_date,
                                0.0625,
                                FrequencyTypes.SEMI_ANNUAL,
                                DayCountTypes.ACT_365F)

num_time_steps = 200


def test_bk_european_exercise():
    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.EUROPEAN

    bermudan_swaption_pay = IborBermudanSwaption(settlement_date,
                                                 exercise_date,
                                                 swapMaturityDate,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swapFixedCoupon,
                                                 swapFixedFrequencyType,
                                                 swapFixedDayCountType)

    fixed_leg_type = SwapTypes.RECEIVE
    exercise_type = FinExerciseTypes.EUROPEAN

    bermudan_swaption_rec = IborBermudanSwaption(settlement_date,
                                                 exercise_date,
                                                 swapMaturityDate,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swapFixedCoupon,
                                                 swapFixedFrequencyType,
                                                 swapFixedDayCountType)

    # Used BK with constant short-rate volatility
    sigma = 0.000001
    a = 0.01
    model = BKTree(sigma, a, num_time_steps)

    valuePay = bermudan_swaption_pay.value(valuation_date, libor_curve, model)
    assert round(valuePay, 4) == 6313.7455

    valueRec = bermudan_swaption_rec.value(valuation_date, libor_curve, model)
    assert valueRec == 0.0

    # Used BK with constant short-rate volatility
    sigma = 0.2
    a = 0.01
    model = BKTree(sigma, a, num_time_steps)

    valuePay = bermudan_swaption_pay.value(valuation_date, libor_curve, model)
    assert round(valuePay, 4) == 15706.6985

    valueRec = bermudan_swaption_rec.value(valuation_date, libor_curve, model)
    assert round(valueRec, 4) == 9392.9531


def test_bk_bermudan_exercise():
    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.BERMUDAN

    bermudan_swaption_pay = IborBermudanSwaption(settlement_date,
                                                 exercise_date,
                                                 swapMaturityDate,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swapFixedCoupon,
                                                 swapFixedFrequencyType,
                                                 swapFixedDayCountType)

    fixed_leg_type = SwapTypes.RECEIVE
    exercise_type = FinExerciseTypes.BERMUDAN

    bermudan_swaption_rec = IborBermudanSwaption(settlement_date,
                                                 exercise_date,
                                                 swapMaturityDate,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swapFixedCoupon,
                                                 swapFixedFrequencyType,
                                                 swapFixedDayCountType)

    # Used BK with constant short-rate volatility
    sigma = 0.000001
    a = 0.01
    model = BKTree(sigma, a, num_time_steps)

    valuePay = bermudan_swaption_pay.value(valuation_date, libor_curve, model)
    assert round(valuePay, 4) == 6313.7455

    valueRec = bermudan_swaption_rec.value(valuation_date, libor_curve, model)
    assert valueRec == 0.0

    # Used BK with constant short-rate volatility
    sigma = 0.20
    a = 0.01
    model = BKTree(sigma, a, num_time_steps)

    valuePay = bermudan_swaption_pay.value(valuation_date, libor_curve, model)
    assert round(valuePay, 4) == 19175.5406

    valueRec = bermudan_swaption_rec.value(valuation_date, libor_curve, model)
    assert round(valueRec, 4) == 12956.6057


def test_bdt_european_exercise():
    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.EUROPEAN

    bermudan_swaption_pay = IborBermudanSwaption(settlement_date,
                                                 exercise_date,
                                                 swapMaturityDate,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swapFixedCoupon,
                                                 swapFixedFrequencyType,
                                                 swapFixedDayCountType)

    fixed_leg_type = SwapTypes.RECEIVE
    exercise_type = FinExerciseTypes.EUROPEAN

    bermudan_swaption_rec = IborBermudanSwaption(settlement_date,
                                                 exercise_date,
                                                 swapMaturityDate,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swapFixedCoupon,
                                                 swapFixedFrequencyType,
                                                 swapFixedDayCountType)

    sigma = 0.00001
    model = BDTTree(sigma, num_time_steps)

    valuePay = bermudan_swaption_pay.value(valuation_date, libor_curve, model)
    assert round(valuePay, 4) == 6313.5155

    valueRec = bermudan_swaption_rec.value(valuation_date, libor_curve, model)
    assert valueRec == 0.0

    sigma = 0.20
    model = BDTTree(sigma, num_time_steps)

    valuePay = bermudan_swaption_pay.value(valuation_date, libor_curve, model)
    assert round(valuePay, 4) == 15969.8968

    valueRec = bermudan_swaption_rec.value(valuation_date, libor_curve, model)
    assert round(valueRec, 4) == 9652.6912


def test_bdt_bermudan_exercise():
    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.BERMUDAN

    bermudan_swaption_pay = IborBermudanSwaption(settlement_date,
                                                 exercise_date,
                                                 swapMaturityDate,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swapFixedCoupon,
                                                 swapFixedFrequencyType,
                                                 swapFixedDayCountType)

    fixed_leg_type = SwapTypes.RECEIVE
    exercise_type = FinExerciseTypes.BERMUDAN

    bermudan_swaption_rec = IborBermudanSwaption(settlement_date,
                                                 exercise_date,
                                                 swapMaturityDate,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swapFixedCoupon,
                                                 swapFixedFrequencyType,
                                                 swapFixedDayCountType)

    sigma = 0.00001
    model = BDTTree(sigma, num_time_steps)

    valuePay = bermudan_swaption_pay.value(valuation_date, libor_curve, model)
    assert round(valuePay, 4) == 6313.5155

    valueRec = bermudan_swaption_rec.value(valuation_date, libor_curve, model)
    assert valueRec == 0.0

    sigma = 0.20
    model = BDTTree(sigma, num_time_steps)

    valuePay = bermudan_swaption_pay.value(valuation_date, libor_curve, model)
    assert round(valuePay, 4) == 19446.0956

    valueRec = bermudan_swaption_rec.value(valuation_date, libor_curve, model)
    assert round(valueRec, 4) == 13256.8421


def test_hw_european_exercise():
    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.EUROPEAN

    bermudan_swaption_pay = IborBermudanSwaption(settlement_date,
                                                 exercise_date,
                                                 swapMaturityDate,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swapFixedCoupon,
                                                 swapFixedFrequencyType,
                                                 swapFixedDayCountType)

    fixed_leg_type = SwapTypes.RECEIVE
    exercise_type = FinExerciseTypes.EUROPEAN

    bermudan_swaption_rec = IborBermudanSwaption(settlement_date,
                                                 exercise_date,
                                                 swapMaturityDate,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swapFixedCoupon,
                                                 swapFixedFrequencyType,
                                                 swapFixedDayCountType)

    sigma = 0.000001
    a = 0.01
    model = HWTree(sigma, a, num_time_steps)

    valuePay = bermudan_swaption_pay.value(valuation_date, libor_curve, model)
    assert round(valuePay, 4) == 6353.3815

    valueRec = bermudan_swaption_rec.value(valuation_date, libor_curve, model)
    assert valueRec == 0.0

    sigma = 0.01
    a = 0.01
    model = HWTree(sigma, a, num_time_steps)

    valuePay = bermudan_swaption_pay.value(valuation_date, libor_curve, model)
    assert round(valuePay, 4) == 13698.0155

    valueRec = bermudan_swaption_rec.value(valuation_date, libor_curve, model)
    assert round(valueRec, 4) == 7344.6339


def test_hw_bermudan_exercise():
    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.BERMUDAN

    bermudan_swaption_pay = IborBermudanSwaption(settlement_date,
                                                 exercise_date,
                                                 swapMaturityDate,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swapFixedCoupon,
                                                 swapFixedFrequencyType,
                                                 swapFixedDayCountType)

    fixed_leg_type = SwapTypes.RECEIVE
    exercise_type = FinExerciseTypes.BERMUDAN

    bermudan_swaption_rec = IborBermudanSwaption(settlement_date,
                                                 exercise_date,
                                                 swapMaturityDate,
                                                 fixed_leg_type,
                                                 exercise_type,
                                                 swapFixedCoupon,
                                                 swapFixedFrequencyType,
                                                 swapFixedDayCountType)

    sigma = 0.000001
    a = 0.01
    model = HWTree(sigma, a, num_time_steps)

    valuePay = bermudan_swaption_pay.value(valuation_date, libor_curve, model)
    assert round(valuePay, 4) == 6353.3815

    valueRec = bermudan_swaption_rec.value(valuation_date, libor_curve, model)
    assert valueRec == 0.0

    sigma = 0.01
    a = 0.01
    model = HWTree(sigma, a, num_time_steps)

    valuePay = bermudan_swaption_pay.value(valuation_date, libor_curve, model)
    assert round(valuePay, 4) == 16609.3646

    valueRec = bermudan_swaption_rec.value(valuation_date, libor_curve, model)
    assert round(valueRec, 4) == 10406.4558
