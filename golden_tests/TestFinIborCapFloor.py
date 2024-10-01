###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

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
from financepy.utils.global_vars import g_days_in_year
from financepy.market.volatility.ibor_cap_vol_curve import IborCapVolCurve
from FinTestCases import FinTestCases, globalTestCaseMode


test_cases = FinTestCases(__file__, globalTestCaseMode)

##############################################################################


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


def test_IborCapFloor():

    todayDate = Date(20, 6, 2019)
    value_dt = todayDate
    start_dt = todayDate.add_weekdays(2)
    maturity_dt = start_dt.add_tenor("1Y")
    libor_curve = test_ibor_depositsAndSwaps(todayDate)

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
        capFloorType = FinCapFloorTypes.CAP
        capfloor = IborCapFloor(start_dt, maturity_dt, capFloorType, k)
        cvalue1 = capfloor.value(value_dt, libor_curve, model1)
        cvalue2 = capfloor.value(value_dt, libor_curve, model2)
        cvalue3 = capfloor.value(value_dt, libor_curve, model3)
        cvalue4 = capfloor.value(value_dt, libor_curve, model4)
        cvalue5 = capfloor.value(value_dt, libor_curve, model5)
        cvalue6 = capfloor.value(value_dt, libor_curve, model6)
        test_cases.print(
            "CAP", k, cvalue1, cvalue2, cvalue3, cvalue4, cvalue5, cvalue6
        )

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
        capFloorType = FinCapFloorTypes.FLOOR
        capfloor = IborCapFloor(start_dt, maturity_dt, capFloorType, k)
        fvalue1 = capfloor.value(value_dt, libor_curve, model1)
        fvalue2 = capfloor.value(value_dt, libor_curve, model2)
        fvalue3 = capfloor.value(value_dt, libor_curve, model3)
        fvalue4 = capfloor.value(value_dt, libor_curve, model4)
        fvalue5 = capfloor.value(value_dt, libor_curve, model5)
        fvalue6 = capfloor.value(value_dt, libor_curve, model6)
        test_cases.print(
            "FLR", k, fvalue1, fvalue2, fvalue3, fvalue4, fvalue5, fvalue6
        )

    ###############################################################################
    # PUT CALL CHECK
    ###############################################################################

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
        capFloorType = FinCapFloorTypes.CAP
        capfloor = IborCapFloor(start_dt, maturity_dt, capFloorType, k)
        cvalue1 = capfloor.value(value_dt, libor_curve, model1)
        cvalue2 = capfloor.value(value_dt, libor_curve, model2)
        cvalue3 = capfloor.value(value_dt, libor_curve, model3)
        cvalue4 = capfloor.value(value_dt, libor_curve, model4)
        cvalue5 = capfloor.value(value_dt, libor_curve, model5)
        cvalue6 = capfloor.value(value_dt, libor_curve, model6)

        capFloorType = FinCapFloorTypes.FLOOR
        capfloor = IborCapFloor(start_dt, maturity_dt, capFloorType, k)
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


###############################################################################


def test_IborCapFloorVolCurve():
    """Aim here is to price cap and caplets using cap and caplet vols and to
    demonstrate they are the same - NOT SURE THAT HULLS BOOKS FORMULA WORKS FOR
    OPTIONS."""

    todayDate = Date(20, 6, 2019)
    value_dt = todayDate
    maturity_dt = value_dt.add_tenor("3Y")
    dc_type = DayCountTypes.THIRTY_E_360
    frequency = FrequencyTypes.ANNUAL

    k = 0.04
    capFloorType = FinCapFloorTypes.CAP
    capFloor = IborCapFloor(
        value_dt, maturity_dt, capFloorType, k, None, frequency, dc_type
    )

    capVolDates = Schedule(
        value_dt, value_dt.add_tenor("10Y"), frequency
    ).generate()

    flat_rate = 0.04
    libor_curve = DiscountCurveFlat(value_dt, flat_rate, frequency, dc_type)

    flat = False
    if flat is True:
        capVolatilities = [20.0] * 11
        capVolatilities[0] = 0.0
    else:
        capVolatilities = [
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

    capVolatilities = np.array(capVolatilities) / 100.0
    capVolatilities[0] = 0.0

    vol_curve = IborCapVolCurve(
        value_dt, capVolDates, capVolatilities, dc_type
    )

    #    print(vol_curve._capletGammas)

    # Value cap using a single flat cap volatility
    tcap = (maturity_dt - value_dt) / g_days_in_year
    vol = vol_curve.cap_vol(maturity_dt)
    model = Black(vol)
    valueCap = capFloor.value(value_dt, libor_curve, model)
    #    print("CAP T", tcap, "VOL:", vol, "VALUE OF CAP:", valueCap)

    # Value cap by breaking it down into caplets using caplet vols
    vCaplets = 0.0
    capletstart_dt = capFloor.capFloorLetDates[1]
    test_cases.header("START", "END", "VOL", "VALUE")

    for caplet_end_dt in capFloor.capFloorLetDates[2:]:
        vol = vol_curve.caplet_vol(caplet_end_dt)
        modelCaplet = Black(vol)
        vCaplet = capFloor.value_caplet_floor_let(
            value_dt, capletstart_dt, caplet_end_dt, libor_curve, modelCaplet
        )

        vCaplets += vCaplet
        test_cases.print(
            "%12s" % capletstart_dt,
            "%s" % caplet_end_dt,
            "%9.5f" % (vol * 100.0),
            "%9.5f" % vCaplet,
        )

        capletstart_dt = caplet_end_dt

    test_cases.header("LABEL", "VALUE")
    test_cases.print("CAPLETS->CAP: ", vCaplets)


###############################################################################


def test_IborCapletHull():

    #  Hull Page 703, example 29.3
    todayDate = Date(20, 6, 2019)
    value_dt = todayDate
    maturity_dt = value_dt.add_tenor("2Y")
    libor_curve = DiscountCurveFlat(
        value_dt, 0.070, FrequencyTypes.QUARTERLY, DayCountTypes.THIRTY_E_360
    )

    k = 0.08
    capFloorType = FinCapFloorTypes.CAP
    capFloor = IborCapFloor(
        value_dt,
        maturity_dt,
        capFloorType,
        k,
        None,
        FrequencyTypes.QUARTERLY,
        DayCountTypes.THIRTY_E_360,
    )

    # Value cap using a single flat cap volatility
    model = Black(0.20)
    capFloor.value(value_dt, libor_curve, model)

    # Value cap by breaking it down into caplets using caplet vols
    capletstart_dt = value_dt.add_tenor("1Y")
    caplet_end_dt = capletstart_dt.add_tenor("3M")

    vCaplet = capFloor.value_caplet_floor_let(
        value_dt, capletstart_dt, caplet_end_dt, libor_curve, model
    )

    # Cannot match Hull due to dates being adjusted
    test_cases.header("CORRECT PRICE", "MODEL_PRICE")
    test_cases.print(517.29, vCaplet)


###############################################################################


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
    option_type = FinCapFloorTypes.CAP
    strike_rate = 0.02

    cap = IborCapFloor(
        start_dt,
        end_dt,
        option_type,
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
    numRepeats = 10
    for i in range(0, numRepeats):
        v = cap.value(value_dt, discount_curve, model)

    end = time.time()
    period = end - start


#    print(v, period/numRepeats)

###############################################################################


test_IborCapletHull()
test_IborCapFloorVolCurve()
test_IborCapFloor()
test_IborCapFloorQLExample()
test_cases.compareTestCases()
