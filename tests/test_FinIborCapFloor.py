###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

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
from financepy.utils.global_types import FinCapFloorTypes


def build_curve(valuation_date):

    depoBasis = DayCountTypes.THIRTY_E_360_ISDA
    depos = []

    spot_days = 0
    settlement_date = valuation_date.add_weekdays(spot_days)
    deposit_rate = 0.05

    depo1 = IborDeposit(settlement_date, "1M", deposit_rate, depoBasis)
    depo2 = IborDeposit(settlement_date, "3M", deposit_rate, depoBasis)
    depo3 = IborDeposit(settlement_date, "6M", deposit_rate, depoBasis)

    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)

    fras = []

    swaps = []
    fixedBasis = DayCountTypes.ACT_365F
    fixedFreq = FrequencyTypes.SEMI_ANNUAL
    fixed_leg_type = SwapTypes.PAY

    swap_rate = 0.05
    swap1 = IborSwap(settlement_date, "1Y", fixed_leg_type,
                     swap_rate, fixedFreq, fixedBasis)
    swap2 = IborSwap(settlement_date, "3Y", fixed_leg_type,
                     swap_rate, fixedFreq, fixedBasis)
    swap3 = IborSwap(settlement_date, "5Y", fixed_leg_type,
                     swap_rate, fixedFreq, fixedBasis)

    swaps.append(swap1)
    swaps.append(swap2)
    swaps.append(swap3)

    libor_curve = IborSingleCurve(valuation_date, depos, fras, swaps)

    return libor_curve


todayDate = Date(20, 6, 2019)
valuation_date = todayDate
start_date = todayDate.add_weekdays(2)
maturity_date = start_date.add_tenor("1Y")
libor_curve = build_curve(todayDate)

model1 = Black(0.20)
model2 = BlackShifted(0.25, 0.0)
model3 = SABR(0.013, 0.5, 0.5, 0.5)
model4 = SABRShifted(0.013, 0.5, 0.5, 0.5, -0.008)
model5 = HWTree(0.30, 0.01)
model6 = Bachelier(0.01)


def test_cap():
    capFloorType = FinCapFloorTypes.CAP

    k = 0.02
    capfloor = IborCapFloor(start_date, maturity_date, capFloorType, k)
    cvalue1 = capfloor.value(valuation_date, libor_curve, model1)
    cvalue2 = capfloor.value(valuation_date, libor_curve, model2)
    cvalue3 = capfloor.value(valuation_date, libor_curve, model3)
    cvalue4 = capfloor.value(valuation_date, libor_curve, model4)
    cvalue5 = capfloor.value(valuation_date, libor_curve, model5)
    cvalue6 = capfloor.value(valuation_date, libor_curve, model6)
    assert round(cvalue1, 4) == 28889.2445
    assert round(cvalue2, 4) == 28889.2482
    assert round(cvalue3, 4) == 28889.2445
    assert round(cvalue4, 4) == 28889.2445
    assert round(cvalue5, 4) == 82406.6040
    assert round(cvalue6, 4) == 28889.3760

    k = 0.05
    capfloor = IborCapFloor(start_date, maturity_date, capFloorType, k)
    cvalue1 = capfloor.value(valuation_date, libor_curve, model1)
    cvalue2 = capfloor.value(valuation_date, libor_curve, model2)
    cvalue3 = capfloor.value(valuation_date, libor_curve, model3)
    cvalue4 = capfloor.value(valuation_date, libor_curve, model4)
    cvalue5 = capfloor.value(valuation_date, libor_curve, model5)
    cvalue6 = capfloor.value(valuation_date, libor_curve, model6)
    assert round(cvalue1, 4) == 1904.9614
    assert round(cvalue2, 4) == 2399.1333
    assert round(cvalue3, 4) == 516.9324
    assert round(cvalue4, 4) == 569.7861
    assert round(cvalue5, 4) == 63709.5424
    assert round(cvalue6, 4) == 1910.0023

    k = 0.08
    capfloor = IborCapFloor(start_date, maturity_date, capFloorType, k)
    cvalue1 = capfloor.value(valuation_date, libor_curve, model1)
    cvalue2 = capfloor.value(valuation_date, libor_curve, model2)
    cvalue3 = capfloor.value(valuation_date, libor_curve, model3)
    cvalue4 = capfloor.value(valuation_date, libor_curve, model4)
    cvalue5 = capfloor.value(valuation_date, libor_curve, model5)
    cvalue6 = capfloor.value(valuation_date, libor_curve, model6)
    assert round(cvalue1, 4) == 3.0923
    assert round(cvalue2, 4) == 21.2585
    assert round(cvalue3, 4) == 0.0023
    assert round(cvalue4, 4) == 0.0187
    assert round(cvalue5, 4) == 53647.5908
    assert round(cvalue6, 4) == 0.1578


def test_floor():
    capFloorType = FinCapFloorTypes.FLOOR

    k = 0.02
    capfloor = IborCapFloor(start_date, maturity_date, capFloorType, k)
    cvalue1 = capfloor.value(valuation_date, libor_curve, model1)
    cvalue2 = capfloor.value(valuation_date, libor_curve, model2)
    cvalue3 = capfloor.value(valuation_date, libor_curve, model3)
    cvalue4 = capfloor.value(valuation_date, libor_curve, model4)
    cvalue5 = capfloor.value(valuation_date, libor_curve, model5)
    cvalue6 = capfloor.value(valuation_date, libor_curve, model6)
    assert round(cvalue1, 4) == 0.0
    assert round(cvalue2, 4) == 0.0037
    assert round(cvalue3, 4) == 0.0
    assert round(cvalue4, 4) == 0.0
    assert round(cvalue5, 4) == 51868.1540
    assert round(cvalue6, 4) == 0.1316

    k = 0.05
    capfloor = IborCapFloor(start_date, maturity_date, capFloorType, k)
    cvalue1 = capfloor.value(valuation_date, libor_curve, model1)
    cvalue2 = capfloor.value(valuation_date, libor_curve, model2)
    cvalue3 = capfloor.value(valuation_date, libor_curve, model3)
    cvalue4 = capfloor.value(valuation_date, libor_curve, model4)
    cvalue5 = capfloor.value(valuation_date, libor_curve, model5)
    cvalue6 = capfloor.value(valuation_date, libor_curve, model6)
    assert round(cvalue1, 4) == 2089.3995
    assert round(cvalue2, 4) == 2583.5715
    assert round(cvalue3, 4) == 701.3705
    assert round(cvalue4, 4) == 754.2243
    assert round(cvalue5, 4) == 62244.0904
    assert round(cvalue6, 4) == 2094.4405

    k = 0.08
    capfloor = IborCapFloor(start_date, maturity_date, capFloorType, k)
    cvalue1 = capfloor.value(valuation_date, libor_curve, model1)
    cvalue2 = capfloor.value(valuation_date, libor_curve, model2)
    cvalue3 = capfloor.value(valuation_date, libor_curve, model3)
    cvalue4 = capfloor.value(valuation_date, libor_curve, model4)
    cvalue5 = capfloor.value(valuation_date, libor_curve, model5)
    cvalue6 = capfloor.value(valuation_date, libor_curve, model6)
    assert round(cvalue1, 4) == 29261.2132
    assert round(cvalue2, 4) == 29279.3794
    assert round(cvalue3, 4) == 29258.1231
    assert round(cvalue4, 4) == 29258.1395
    assert round(cvalue5, 4) == 81255.1368
    assert round(cvalue6, 4) == 29258.2786
