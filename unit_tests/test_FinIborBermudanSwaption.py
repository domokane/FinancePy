# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.bdt_tree import BDTTree
from financepy.models.hw_tree import HWTree
from financepy.models.bk_tree import BKTree
from financepy.products.rates.ibor_bermudan_swaption import (
    IborBermudanSwaption,
)
from financepy.utils.global_types import FinExerciseTypes
from financepy.utils.global_types import SwapTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date


value_dt = Date(1, 1, 2011)
settle_dt = value_dt
exercise_dt = settle_dt.add_years(1)
swap_maturity_dt = settle_dt.add_years(4)

swap_fixed_cpn = 0.060
swap_fixed_freq_type = FrequencyTypes.SEMI_ANNUAL
swap_fixed_day_count_type = DayCountTypes.ACT_365F

libor_curve = DiscountCurveFlat(
    value_dt, 0.0625, FrequencyTypes.SEMI_ANNUAL, DayCountTypes.ACT_365F
)

num_time_steps = 200

########################################################################################


def test_bk_european_exercise():

    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.EUROPEAN

    bermudan_swaption_pay = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    fixed_leg_type = SwapTypes.RECEIVE
    exercise_type = FinExerciseTypes.EUROPEAN

    bermudan_swaption_rec = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    # Used BK with constant short-rate volatility
    sigma = 0.000001
    a = 0.01
    model = BKTree(sigma, a, num_time_steps)

    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    assert round(value_pay, 4) == 6313.7455

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    assert value_rec == 0.0

    # Used BK with constant short-rate volatility
    sigma = 0.2
    a = 0.01
    model = BKTree(sigma, a, num_time_steps)

    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    assert round(value_pay, 4) == 15706.6985

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    assert round(value_rec, 4) == 9392.9531


########################################################################################


def test_bk_bermudan_exercise():

    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.BERMUDAN

    bermudan_swaption_pay = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    fixed_leg_type = SwapTypes.RECEIVE
    exercise_type = FinExerciseTypes.BERMUDAN

    bermudan_swaption_rec = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    # Used BK with constant short-rate volatility
    sigma = 0.000001
    a = 0.01
    model = BKTree(sigma, a, num_time_steps)

    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    assert round(value_pay, 4) == 6313.7455

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    assert value_rec == 0.0

    # Used BK with constant short-rate volatility
    sigma = 0.20
    a = 0.01
    model = BKTree(sigma, a, num_time_steps)

    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    assert round(value_pay, 4) == 19175.5406

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    assert round(value_rec, 4) == 12956.6057


########################################################################################


def test_bdt_european_exercise():

    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.EUROPEAN

    bermudan_swaption_pay = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    fixed_leg_type = SwapTypes.RECEIVE
    exercise_type = FinExerciseTypes.EUROPEAN

    bermudan_swaption_rec = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    sigma = 0.00001
    model = BDTTree(sigma, num_time_steps)

    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    assert round(value_pay, 4) == 6313.5155

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    assert value_rec == 0.0

    sigma = 0.20
    model = BDTTree(sigma, num_time_steps)

    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    assert round(value_pay, 4) == 15969.8968

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    assert round(value_rec, 4) == 9652.6912


########################################################################################


def test_bdt_bermudan_exercise():

    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.BERMUDAN

    bermudan_swaption_pay = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    fixed_leg_type = SwapTypes.RECEIVE
    exercise_type = FinExerciseTypes.BERMUDAN

    bermudan_swaption_rec = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    sigma = 0.00001
    model = BDTTree(sigma, num_time_steps)

    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    assert round(value_pay, 4) == 6313.5155

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    assert value_rec == 0.0

    sigma = 0.20
    model = BDTTree(sigma, num_time_steps)

    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    assert round(value_pay, 4) == 19446.0956

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    assert round(value_rec, 4) == 13256.8421


########################################################################################


def test_hw_european_exercise():

    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.EUROPEAN

    bermudan_swaption_pay = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    fixed_leg_type = SwapTypes.RECEIVE
    exercise_type = FinExerciseTypes.EUROPEAN

    bermudan_swaption_rec = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    sigma = 0.000001
    a = 0.01
    model = HWTree(sigma, a, num_time_steps)

    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    assert round(value_pay, 4) == 6353.3815

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    assert value_rec == 0.0

    sigma = 0.01
    a = 0.01
    model = HWTree(sigma, a, num_time_steps)

    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    assert round(value_pay, 4) == 13698.0155

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    assert round(value_rec, 4) == 7344.6339


########################################################################################


def test_hw_bermudan_exercise():

    fixed_leg_type = SwapTypes.PAY
    exercise_type = FinExerciseTypes.BERMUDAN

    bermudan_swaption_pay = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    fixed_leg_type = SwapTypes.RECEIVE
    exercise_type = FinExerciseTypes.BERMUDAN

    bermudan_swaption_rec = IborBermudanSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        fixed_leg_type,
        exercise_type,
        swap_fixed_cpn,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    sigma = 0.000001
    a = 0.01
    model = HWTree(sigma, a, num_time_steps)

    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    assert round(value_pay, 4) == 6353.3815

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    assert value_rec == 0.0

    sigma = 0.01
    a = 0.01
    model = HWTree(sigma, a, num_time_steps)

    value_pay = bermudan_swaption_pay.value(value_dt, libor_curve, model)
    assert round(value_pay, 4) == 16609.3646

    value_rec = bermudan_swaption_rec.value(value_dt, libor_curve, model)
    assert round(value_rec, 4) == 10406.4558
