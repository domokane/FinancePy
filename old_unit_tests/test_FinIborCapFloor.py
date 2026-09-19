# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.models.hw_tree import HWTree
from financepy.models.sabr_shifted import SABRShifted
from financepy.models.sabr import SABR
from financepy.models.black_shifted import BlackShifted
from financepy.models.bachelier import Bachelier
from financepy.models.black import Black
from financepy.utils.global_types import SwapTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.rates.ibor_cap_floor import IborCapFloor
from financepy.utils.global_types import CapFloorTypes

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


today_date = Date(20, 6, 2019)
value_dt = today_date
start_dt = today_date.add_weekdays(2)
maturity_dt = start_dt.add_tenor("1Y")
libor_curve = build_curve(today_date)

model1 = Black(0.20)
model2 = BlackShifted(0.25, 0.0)
model3 = SABR(0.013, 0.5, 0.5, 0.5)
model4 = SABRShifted(0.013, 0.5, 0.5, 0.5, -0.008)
model5 = HWTree(0.30, 0.01)
model6 = Bachelier(0.01)

########################################################################################


def test_cap():

    cap_floor_type = CapFloorTypes.CAP

    k = 0.02
    capfloor = IborCapFloor(start_dt, maturity_dt, cap_floor_type, k)
    cvalue1 = capfloor.value(value_dt, libor_curve, model1)
    cvalue2 = capfloor.value(value_dt, libor_curve, model2)
    cvalue3 = capfloor.value(value_dt, libor_curve, model3)
    cvalue4 = capfloor.value(value_dt, libor_curve, model4)
    cvalue5 = capfloor.value(value_dt, libor_curve, model5)
    cvalue6 = capfloor.value(value_dt, libor_curve, model6)
    assert round(cvalue1, 4) == 28889.4749
    assert round(cvalue2, 4) == 28889.4786
    assert round(cvalue3, 4) == 28889.4749
    assert round(cvalue4, 4) == 28889.4749
    assert round(cvalue5, 4) == 82372.5602
    assert round(cvalue6, 4) == 28889.6059

    k = 0.05
    capfloor = IborCapFloor(start_dt, maturity_dt, cap_floor_type, k)
    cvalue1 = capfloor.value(value_dt, libor_curve, model1)
    cvalue2 = capfloor.value(value_dt, libor_curve, model2)
    cvalue3 = capfloor.value(value_dt, libor_curve, model3)
    cvalue4 = capfloor.value(value_dt, libor_curve, model4)
    cvalue5 = capfloor.value(value_dt, libor_curve, model5)
    cvalue6 = capfloor.value(value_dt, libor_curve, model6)
    assert round(cvalue1, 0) == 1905
    assert round(cvalue2, 0) == 2399
    assert round(cvalue3, 0) == 517
    assert round(cvalue4, 0) == 570
    assert round(cvalue5, 0) == 63678
    assert round(cvalue6, 0) == 1910

    k = 0.08
    capfloor = IborCapFloor(start_dt, maturity_dt, cap_floor_type, k)
    cvalue1 = capfloor.value(value_dt, libor_curve, model1)
    cvalue2 = capfloor.value(value_dt, libor_curve, model2)
    cvalue3 = capfloor.value(value_dt, libor_curve, model3)
    cvalue4 = capfloor.value(value_dt, libor_curve, model4)
    cvalue5 = capfloor.value(value_dt, libor_curve, model5)
    cvalue6 = capfloor.value(value_dt, libor_curve, model6)


    assert round(cvalue1, 4) == 3.1029
    assert round(cvalue2, 4) == 21.3032
    assert round(cvalue3, 4) == 0.0023
    assert round(cvalue4, 4) == 0.0188
    assert round(cvalue5, 0) == 53619
    assert round(cvalue6, 4) == 0.1585


########################################################################################


def test_floor():

    cap_floor_type = CapFloorTypes.FLOOR

    k = 0.02
    capfloor = IborCapFloor(start_dt, maturity_dt, cap_floor_type, k)
    cvalue1 = capfloor.value(value_dt, libor_curve, model1)
    cvalue2 = capfloor.value(value_dt, libor_curve, model2)
    cvalue3 = capfloor.value(value_dt, libor_curve, model3)
    cvalue4 = capfloor.value(value_dt, libor_curve, model4)
    cvalue5 = capfloor.value(value_dt, libor_curve, model5)
    cvalue6 = capfloor.value(value_dt, libor_curve, model6)
    assert round(cvalue1, 4) == 0.0
    assert round(cvalue2, 4) == 0.0037
    assert round(cvalue3, 4) == 0.0
    assert round(cvalue4, 4) == 0.0
    assert round(cvalue5, 0) == 51899
    assert round(cvalue6, 4) == 0.1310

    k = 0.05
    capfloor = IborCapFloor(start_dt, maturity_dt, cap_floor_type, k)
    cvalue1 = capfloor.value(value_dt, libor_curve, model1)
    cvalue2 = capfloor.value(value_dt, libor_curve, model2)
    cvalue3 = capfloor.value(value_dt, libor_curve, model3)
    cvalue4 = capfloor.value(value_dt, libor_curve, model4)
    cvalue5 = capfloor.value(value_dt, libor_curve, model5)
    cvalue6 = capfloor.value(value_dt, libor_curve, model6)
    assert round(cvalue1, 0) == 2089
    assert round(cvalue2, 0) == 2584
    assert round(cvalue3, 0) == 702
    assert round(cvalue4, 0) == 754
    assert round(cvalue5, 0) == 62278
    assert round(cvalue6, 0) == 2094

    k = 0.08
    capfloor = IborCapFloor(start_dt, maturity_dt, cap_floor_type, k)
    cvalue1 = capfloor.value(value_dt, libor_curve, model1)
    cvalue2 = capfloor.value(value_dt, libor_curve, model2)
    cvalue3 = capfloor.value(value_dt, libor_curve, model3)
    cvalue4 = capfloor.value(value_dt, libor_curve, model4)
    cvalue5 = capfloor.value(value_dt, libor_curve, model5)
    cvalue6 = capfloor.value(value_dt, libor_curve, model6)
    assert round(cvalue1, 0) == 29261
    assert round(cvalue2, 0) == 29279
    assert round(cvalue3, 0) == 29258
    assert round(cvalue4, 0) == 29258
    assert round(cvalue5, 0) == 81293
    assert round(cvalue6, 0) == 29258
