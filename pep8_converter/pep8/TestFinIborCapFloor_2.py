# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
import time

import sys

sys.path.append("..")

from financepy.utils.schedule import Schedule
from financepy.utils.global_types import FinCapFloorTypes
from financepy.products.rates.ibor_cap_floor import IborCapFloor
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
from financepy.utils.calendar import CalendarTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.global_types import SwapTypes
from financepy.market.curves.discount_curve_zeros import DiscountCurveZeros
from financepy.market.curves.interpolator import InterpTypes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black import Black
from financepy.models.bachelier import Bachelier
from financepy.models.black_shifted import BlackShifted
from financepy.models.sabr import SABR
from financepy.models.sabr_shifted import SABRShifted
from financepy.models.hw_tree import HWTree
from financepy.utils.global_vars import G_DAYS_IN_YEARS
from financepy.market.volatility.ibor_cap_vol_curve import IborCapVolCurve
from FinTestCases import FinTestCases, global_test_case_mode


test_cases = FinTestCases(__file__, global_test_case_mode)



########################################################################################


def test_ibor_depositsAndSwaps(value_dt):

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
    swap1 = IborSwap(settle_dt, "1Y", fixed_leg_type, swap_rate, fixed_freq, fixed_basis)
    swap2 = IborSwap(settle_dt, "3Y", fixed_leg_type, swap_rate, fixed_freq, fixed_basis)
    swap3 = IborSwap(settle_dt, "5Y", fixed_leg_type, swap_rate, fixed_freq, fixed_basis)

    swaps.append(swap1)
    swaps.append(swap2)
    swaps.append(swap3)

    libor_curve = IborSingleCurve(value_dt, depos, fras, swaps)

    return libor_curve


##########################################################################


########################################################################################


def test_IborCapFloor():

    today_date = Date(20, 6, 2019)
    value_dt = today_date
    start_dt = today_date.add_weekdays(2)
    maturity_dt = start_dt.add_tenor("1Y")
    libor_curve = test_ibor_depositsAndSwaps(today_date)

    # The capfloor has begun
    # last_fixing = 0.028

    ##########################################################################
    # COMPARISON OF MODELS
    ##########################################################################

    strikes = np.linspace(0.02, 0.08, 5)

    test_cases.header(
        "LABEL",
        "STRIKE",
        "BLK",
        "BLK_SHFTD",
        "SABR",
        "SABR_SHFTD",
        "HW",
        "BACH",
    )

    model1 = Black(0.20)
    model2 = BlackShifted(0.25, 0.0)
    model3 = SABR(0.013, 0.5, 0.5, 0.5)
    model4 = SABRShifted(0.013, 0.5, 0.5, 0.5, -0.008)
    model5 = HWTree(0.30, 0.01)
    model6 = Bachelier(0.01)

    for k in strikes:
        cap_floor_type = FinCapFloorTypes.CAP
        capfloor = IborCapFloor(start_dt, maturity_dt, cap_floor_type, k)
        cvalue1 = capfloor.value(value_dt, libor_curve, model1)
        cvalue2 = capfloor.value(value_dt, libor_curve, model2)
        cvalue3 = capfloor.value(value_dt, libor_curve, model3)
        cvalue4 = capfloor.value(value_dt, libor_curve, model4)
        cvalue5 = capfloor.value(value_dt, libor_curve, model5)
        cvalue6 = capfloor.value(value_dt, libor_curve, model6)
        test_cases.print("CAP", k, cvalue1, cvalue2, cvalue3, cvalue4, cvalue5, cvalue6)

    test_cases.header(
        "LABEL",
        "STRIKE",
        "BLK",
        "BLK_SHFTD",
        "SABR",
        "SABR_SHFTD",
        "HW",
        "BACH",
    )

    for k in strikes:
        cap_floor_type = FinCapFloorTypes.FLOOR
        capfloor = IborCapFloor(start_dt, maturity_dt, cap_floor_type, k)
        fvalue1 = capfloor.value(value_dt, libor_curve, model1)
        fvalue2 = capfloor.value(value_dt, libor_curve, model2)
        fvalue3 = capfloor.value(value_dt, libor_curve, model3)
        fvalue4 = capfloor.value(value_dt, libor_curve, model4)
        fvalue5 = capfloor.value(value_dt, libor_curve, model5)
        fvalue6 = capfloor.value(value_dt, libor_curve, model6)
        test_cases.print("FLR", k, fvalue1, fvalue2, fvalue3, fvalue4, fvalue5, fvalue6)

    ########################################################################################
    # PUT CALL CHECK
    ########################################################################################

    test_cases.header(
        "LABEL",
        "STRIKE",
        "BLK",
        "BLK_SHFTD",
        "SABR",
        "SABR SHFTD",
        "HW",
        "BACH",
    )

    for k in strikes:
        cap_floor_type = FinCapFloorTypes.CAP
        capfloor = IborCapFloor(start_dt, maturity_dt, cap_floor_type, k)
        cvalue1 = capfloor.value(value_dt, libor_curve, model1)
        cvalue2 = capfloor.value(value_dt, libor_curve, model2)
        cvalue3 = capfloor.value(value_dt, libor_curve, model3)
        cvalue4 = capfloor.value(value_dt, libor_curve, model4)
        cvalue5 = capfloor.value(value_dt, libor_curve, model5)
        cvalue6 = capfloor.value(value_dt, libor_curve, model6)

        cap_floor_type = FinCapFloorTypes.FLOOR
        capfloor = IborCapFloor(start_dt, maturity_dt, cap_floor_type, k)
        fvalue1 = capfloor.value(value_dt, libor_curve, model1)
        fvalue2 = capfloor.value(value_dt, libor_curve, model2)
        fvalue3 = capfloor.value(value_dt, libor_curve, model3)
        fvalue4 = capfloor.value(value_dt, libor_curve, model4)
        fvalue5 = capfloor.value(value_dt, libor_curve, model5)
        fvalue6 = capfloor.value(value_dt, libor_curve, model6)

        pcvalue1 = cvalue1 - fvalue1
        pcvalue2 = cvalue2 - fvalue2
        pcvalue3 = cvalue3 - fvalue3
        pcvalue4 = cvalue4 - fvalue4
        pcvalue5 = cvalue5 - fvalue5
        pcvalue6 = cvalue6 - fvalue6

        test_cases.print(
            "PUT_CALL",
            k,
            pcvalue1,
            pcvalue2,
            pcvalue3,
            pcvalue4,
            pcvalue5,
            pcvalue6,
        )


########################################################################################


########################################################################################


def test_IborCapFloorVolCurve():
    """aim here is to price cap and caplets using cap and caplet vols and to
    demonstrate they are the same - NOT SURE THAT HULLS BOOKS FORMULA WORKS FOR
    OPTIONS."""

    today_date = Date(20, 6, 2019)
    value_dt = today_date
    maturity_dt = value_dt.add_tenor("3Y")
    dc_type = DayCountTypes.THIRTY_E_360
    frequency = FrequencyTypes.ANNUAL

    k = 0.04
    cap_floor_type = FinCapFloorTypes.CAP
    cap_floor = IborCapFloor(
        value_dt, maturity_dt, cap_floor_type, k, None, frequency, dc_type
    )

    cap_vol_dates = Schedule(value_dt, value_dt.add_tenor("10Y"), frequency).generate()

    flat_rate = 0.04
    libor_curve = DiscountCurveFlat(value_dt, flat_rate, frequency, dc_type)

    flat = False
    if flat is True:
        cap_volatilities = [20.0] * 11
        cap_volatilities[0] = 0.0
    else:
        cap_volatilities = [
            0.00,
            15.50,
            18.25,
            17.91,
            17.74,
            17.27,
            16.79,
            16.30,
            16.01,
            15.76,
            15.54,
        ]

    cap_volatilities = np.array(cap_volatilities) / 100.0
    cap_volatilities[0] = 0.0

    vol_curve = IborCapVolCurve(value_dt, cap_vol_dates, cap_volatilities, dc_type)

    #    print(vol_curve._capletGammas)

    # Value cap using a single flat cap volatility
    tcap = (maturity_dt - value_dt) / G_DAYS_IN_YEARS
    vol = vol_curve.cap_vol(maturity_dt)
    model = Black(vol)
    value_cap = cap_floor.value(value_dt, libor_curve, model)
    #    print("CAP T", tcap, "VOL:", vol, "VALUE OF CAP:", value_cap)

    # Value cap by breaking it down into caplets using caplet vols
    v_caplets = 0.0
    capletstart_dt = cap_floor.caplet_floorlet_dates[1]
    test_cases.header("START", "END", "VOL", "VALUE")

    for caplet_end_dt in cap_floor.caplet_floorlet_dates[2:]:
        vol = vol_curve.caplet_vol(caplet_end_dt)
        model_caplet = Black(vol)
        v_caplet = cap_floor.value_caplet_floor_let(
            value_dt, capletstart_dt, caplet_end_dt, libor_curve, model_caplet
        )

        v_caplets += v_caplet
        test_cases.print(
            "%12s" % capletstart_dt,
            "%s" % caplet_end_dt,
            "%9.5f" % (vol * 100.0),
            "%9.5f" % v_caplet,
        )

        capletstart_dt = caplet_end_dt

    test_cases.header("LABEL", "VALUE")
    test_cases.print("CAPLETS->CAP: ", v_caplets)


########################################################################################


########################################################################################


def test_IborCapletHull():

    #  Hull Page 703, example 29.3
    today_date = Date(20, 6, 2019)
    value_dt = today_date
    maturity_dt = value_dt.add_tenor("2Y")
    libor_curve = DiscountCurveFlat(
        value_dt, 0.070, FrequencyTypes.QUARTERLY, DayCountTypes.THIRTY_E_360
    )

    k = 0.08
    cap_floor_type = FinCapFloorTypes.CAP
    cap_floor = IborCapFloor(
        value_dt,
        maturity_dt,
        cap_floor_type,
        k,
        None,
        FrequencyTypes.QUARTERLY,
        DayCountTypes.THIRTY_E_360,
    )

    # Value cap using a single flat cap volatility
    model = Black(0.20)
    cap_floor.value(value_dt, libor_curve, model)

    # Value cap by breaking it down into caplets using caplet vols
    capletstart_dt = value_dt.add_tenor("1Y")
    caplet_end_dt = capletstart_dt.add_tenor("3M")

    v_caplet = cap_floor.value_caplet_floor_let(
        value_dt, capletstart_dt, caplet_end_dt, libor_curve, model
    )

    # Cannot match Hull due to dates being adjusted
    test_cases.header("CORRECT PRICE", "MODEL_PRICE")
    test_cases.print(517.29, v_caplet)


########################################################################################


########################################################################################


def test_IborCapFloorQLExample():

    value_dt = Date(14, 6, 2016)

    dates = [
        Date(14, 6, 2016),
        Date(14, 9, 2016),
        Date(14, 12, 2016),
        Date(14, 6, 2017),
        Date(14, 6, 2019),
        Date(14, 6, 2021),
        Date(15, 6, 2026),
        Date(16, 6, 2031),
        Date(16, 6, 2036),
        Date(14, 6, 2046),
    ]

    rates = [
        0.000000,
        0.006616,
        0.007049,
        0.007795,
        0.009599,
        0.011203,
        0.015068,
        0.017583,
        0.018998,
        0.020080,
    ]

    freq_type = FrequencyTypes.ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ISDA

    discount_curve = DiscountCurveZeros(
        value_dt,
        dates,
        rates,
        freq_type,
        dc_type,
        InterpTypes.LINEAR_ZERO_RATES,
    )

    start_dt = Date(14, 6, 2016)
    end_dt = Date(14, 6, 2026)
    cal_type = CalendarTypes.UNITED_STATES
    bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    freq_type = FrequencyTypes.QUARTERLY
    dg_type = DateGenRuleTypes.FORWARD
    last_fixing = 0.0065560
    notional = 1000000
    dc_type = DayCountTypes.ACT_360
    opt_type = FinCapFloorTypes.CAP
    strike_rate = 0.02

    cap = IborCapFloor(
        start_dt,
        end_dt,
        opt_type,
        strike_rate,
        last_fixing,
        freq_type,
        dc_type,
        notional,
        cal_type,
        bd_type,
        dg_type,
    )

    black_vol = 0.547295
    model = Black(black_vol)

    start = time.time()
    num_repeats = 10
    for i in range(0, num_repeats):
        v = cap.value(value_dt, discount_curve, model)

    end = time.time()
    period = end - start


#    print(v, period/num_repeats)

########################################################################################


test_IborCapletHull()
test_IborCapFloorVolCurve()
test_IborCapFloor()
test_IborCapFloorQLExample()
test_cases.compare_test_cases()

########################################################################################

