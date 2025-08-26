# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.models.bk_tree import BKTree
from financepy.models.hw_tree import HWTree
from financepy.models.sabr_shifted import SABRShifted
from financepy.models.sabr import SABR
from financepy.models.black_shifted import BlackShifted
from financepy.models.black import Black
from financepy.products.rates.ibor_swaption import SwapTypes
from financepy.products.rates.ibor_swaption import IborSwaption
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date

########################################################################################


def build_curve(value_dt):

    depo_basis = DayCountTypes.THIRTY_E_360_ISDA
    depos = []

    spot_days = 0
    settle_dt = value_dt.add_weekdays(spot_days)
    deposit_rate = 0.05

    depo1 = IborDeposit(settle_dt, "1M", deposit_rate, depo_basis)
    depo2 = IborDeposit(settle_dt, "3M", deposit_rate, depo_basis)
    depo3 = IborDeposit(settle_dt, "6M", deposit_rate, depo_basis)

    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)

    fras = []

    swaps = []
    fixed_basis = DayCountTypes.ACT_365F
    fixed_freq = FrequencyTypes.SEMI_ANNUAL
    fixed_leg_type = SwapTypes.PAY

    swap_rate = 0.05
    swap1 = IborSwap(
        settle_dt, "1Y", fixed_leg_type, swap_rate, fixed_freq, fixed_basis
    )
    swap2 = IborSwap(
        settle_dt, "3Y", fixed_leg_type, swap_rate, fixed_freq, fixed_basis
    )
    swap3 = IborSwap(
        settle_dt, "5Y", fixed_leg_type, swap_rate, fixed_freq, fixed_basis
    )

    swaps.append(swap1)
    swaps.append(swap2)
    swaps.append(swap3)

    libor_curve = IborSingleCurve(value_dt, depos, fras, swaps)

    return libor_curve


value_dt = Date(1, 1, 2011)


exercise_dt = Date(1, 1, 2012)
swap_maturity_dt = Date(1, 1, 2017)

swap_fixed_freq_type = FrequencyTypes.SEMI_ANNUAL
swap_fixed_day_count_type = DayCountTypes.ACT_365F

model1 = Black(0.00001)
model2 = BlackShifted(0.00001, 0.0)
model3 = SABR(0.013, 0.5, 0.5, 0.5)
model4 = SABRShifted(0.013, 0.5, 0.5, 0.5, -0.008)
model5 = HWTree(0.00001, 0.00001)
model6 = BKTree(0.01, 0.01)

settle_dt = value_dt.add_weekdays(2)

libor_curve = build_curve(value_dt)

########################################################################################


def test_pay():

    libor_curve = build_curve(value_dt)
    swaption_type = SwapTypes.PAY

    k = 0.02
    swaption = IborSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        swaption_type,
        k,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    swap1 = swaption.value(value_dt, libor_curve, model1)
    swap2 = swaption.value(value_dt, libor_curve, model2)
    swap3 = swaption.value(value_dt, libor_curve, model3)
    swap4 = swaption.value(value_dt, libor_curve, model4)
    swap5 = swaption.value(value_dt, libor_curve, model5)
    swap6 = swaption.value(value_dt, libor_curve, model6)
    assert round(swap1, 0) == 125087
    assert round(swap2, 0) == 125087
    assert round(swap3, 0) == 125087
    assert round(swap4, 0) == 125087
    assert round(swap5, 0) == 125684
    assert round(swap6, 0) == 124501

    k = 0.035
    swaption = IborSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        swaption_type,
        k,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    swap1 = swaption.value(value_dt, libor_curve, model1)
    swap2 = swaption.value(value_dt, libor_curve, model2)
    swap3 = swaption.value(value_dt, libor_curve, model3)
    swap4 = swaption.value(value_dt, libor_curve, model4)
    swap5 = swaption.value(value_dt, libor_curve, model5)
    swap6 = swaption.value(value_dt, libor_curve, model6)
    assert round(swap1, 1) == 62492.6
    assert round(swap2, 1) == 62492.6
    assert round(swap3, 1) == 62492.6
    assert round(swap4, 1) == 62492.8
    assert round(swap5, 1) == 63098.5
    assert round(swap6, 1) == 62307.2

    k = 0.065
    swaption = IborSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        swaption_type,
        k,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    swap1 = swaption.value(value_dt, libor_curve, model1)
    swap2 = swaption.value(value_dt, libor_curve, model2)
    swap3 = swaption.value(value_dt, libor_curve, model3)
    swap4 = swaption.value(value_dt, libor_curve, model4)
    swap5 = swaption.value(value_dt, libor_curve, model5)
    swap6 = swaption.value(value_dt, libor_curve, model6)
    assert round(swap1, 4) == 0.0
    assert round(swap2, 4) == 0.0
    assert round(swap3, 1) == 22.1
    assert round(swap4, 1) == 60.3
    assert round(swap5, 4) == 0.0
    assert round(swap6, 4) == 0.0


########################################################################################


def test_receive():

    swaption_type = SwapTypes.RECEIVE

    k = 0.02
    swaption = IborSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        swaption_type,
        k,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    swap1 = swaption.value(value_dt, libor_curve, model1)
    swap2 = swaption.value(value_dt, libor_curve, model2)
    swap3 = swaption.value(value_dt, libor_curve, model3)
    swap4 = swaption.value(value_dt, libor_curve, model4)
    swap5 = swaption.value(value_dt, libor_curve, model5)
    swap6 = swaption.value(value_dt, libor_curve, model6)
    assert round(swap1, 4) == 0.0
    assert round(swap2, 4) == 0.0
    assert round(swap3, 4) == 0.0
    assert round(swap4, 4) == 0.0046
    assert round(swap5, 4) == 0.0
    assert round(swap6, 4) == 0.0

    k = 0.05
    swaption = IborSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        swaption_type,
        k,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    swap1 = swaption.value(value_dt, libor_curve, model1)
    swap2 = swaption.value(value_dt, libor_curve, model2)
    swap3 = swaption.value(value_dt, libor_curve, model3)
    swap4 = swaption.value(value_dt, libor_curve, model4)
    swap5 = swaption.value(value_dt, libor_curve, model5)
    swap6 = swaption.value(value_dt, libor_curve, model6)
    assert round(swap1, 1) == 101.8
    assert round(swap2, 1) == 101.8
    assert round(swap3, 1) == 4945.4
    assert round(swap4, 1) == 5392.6
    assert round(swap5, 4) == 0.0
    assert round(swap6, 1) == 762.5

    k = 0.08
    swaption = IborSwaption(
        settle_dt,
        exercise_dt,
        swap_maturity_dt,
        swaption_type,
        k,
        swap_fixed_freq_type,
        swap_fixed_day_count_type,
    )

    swap1 = swaption.value(value_dt, libor_curve, model1)
    swap2 = swaption.value(value_dt, libor_curve, model2)
    swap3 = swaption.value(value_dt, libor_curve, model3)
    swap4 = swaption.value(value_dt, libor_curve, model4)
    swap5 = swaption.value(value_dt, libor_curve, model5)
    swap6 = swaption.value(value_dt, libor_curve, model6)
    assert round(swap1, 1) == 125290.5
    assert round(swap2, 1) == 125290.5
    assert round(swap3, 1) == 125291.1
    assert round(swap4, 1) == 125293.6
    assert round(swap5, 1) == 124657.1
    assert round(swap6, 1) == 124274.9
