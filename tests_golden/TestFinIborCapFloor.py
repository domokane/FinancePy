###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################


import time
import numpy as np

import sys
sys.path.append("..")

from financepy.utils.global_types import FinCapFloorTypes
from financepy.products.rates.ibor_cap_floor import IborCapFloor
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.rates.ibor_swap import FinSwapTypes
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_single_curve import IborSingleCurve

from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date

from financepy.utils.calendar import CalendarTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.calendar import DateGenRuleTypes

from financepy.utils.global_types import FinSwapTypes

from financepy.market.discount.curve_zeros import DiscountCurveZeros
from financepy.market.discount.interpolator import InterpTypes
from financepy.market.discount.curve_flat import DiscountCurveFlat

from financepy.models.black import FinModelBlack
from financepy.models.bachelier import FinModelBachelier
from financepy.models.black_shifted import FinModelBlackShifted
from financepy.models.sabr import FinModelSABR
from financepy.models.sabr_shifted import FinModelSABRShifted
from financepy.models.rates_hull_white_tree import FinModelRatesHW

from financepy.utils.global_vars import gDaysInYear

from financepy.market.volatility.ibor_cap_vol_curve import IborCapVolCurve
from financepy.utils.schedule import Schedule

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##############################################################################

def test_FinIborDepositsAndSwaps(valuation_date):

    depoBasis = DayCountTypes.THIRTY_E_360_ISDA
    depos = []

    spot_days = 0
    settlement_date = valuation_date.addWeekDays(spot_days)
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
    fixed_leg_type = FinSwapTypes.PAY

    swap_rate = 0.05
    swap1 = IborSwap(settlement_date, "1Y", fixed_leg_type, swap_rate, fixedFreq, fixedBasis)
    swap2 = IborSwap(settlement_date, "3Y", fixed_leg_type, swap_rate, fixedFreq, fixedBasis)
    swap3 = IborSwap(settlement_date, "5Y", fixed_leg_type, swap_rate, fixedFreq, fixedBasis)

    swaps.append(swap1)
    swaps.append(swap2)
    swaps.append(swap3)

    libor_curve = IborSingleCurve(valuation_date, depos, fras, swaps)

    return libor_curve

##########################################################################


def test_FinIborCapFloor():

    todayDate = Date(20, 6, 2019)
    valuation_date = todayDate
    start_date = todayDate.addWeekDays(2)
    maturity_date = start_date.addTenor("1Y")
    libor_curve = test_FinIborDepositsAndSwaps(todayDate)

    # The capfloor has begun
    # lastFixing = 0.028

    ##########################################################################
    # COMPARISON OF MODELS
    ##########################################################################

    strikes = np.linspace(0.02, 0.08, 5)

    testCases.header("LABEL", "STRIKE", "BLK", "BLK_SHFTD", "SABR",
                     "SABR_SHFTD", "HW", "BACH")

    model1 = FinModelBlack(0.20)
    model2 = FinModelBlackShifted(0.25, 0.0)
    model3 = FinModelSABR(0.013, 0.5, 0.5, 0.5)
    model4 = FinModelSABRShifted(0.013, 0.5, 0.5, 0.5, -0.008)
    model5 = FinModelRatesHW(0.30, 0.01)
    model6 = FinModelBachelier(0.01)

    for k in strikes:
        capFloorType = FinCapFloorTypes.CAP
        capfloor = IborCapFloor(start_date, maturity_date, capFloorType, k)
        cvalue1 = capfloor.value(valuation_date, libor_curve, model1)
        cvalue2 = capfloor.value(valuation_date, libor_curve, model2)
        cvalue3 = capfloor.value(valuation_date, libor_curve, model3)
        cvalue4 = capfloor.value(valuation_date, libor_curve, model4)
        cvalue5 = capfloor.value(valuation_date, libor_curve, model5)
        cvalue6 = capfloor.value(valuation_date, libor_curve, model6)
        testCases.print("CAP", k, cvalue1, cvalue2, cvalue3, cvalue4, cvalue5, cvalue6)

    testCases.header("LABEL", "STRIKE", "BLK", "BLK_SHFTD", "SABR",
                     "SABR_SHFTD", "HW", "BACH")

    for k in strikes:
        capFloorType = FinCapFloorTypes.FLOOR
        capfloor = IborCapFloor(start_date, maturity_date, capFloorType, k)
        fvalue1 = capfloor.value(valuation_date, libor_curve, model1)
        fvalue2 = capfloor.value(valuation_date, libor_curve, model2)
        fvalue3 = capfloor.value(valuation_date, libor_curve, model3)
        fvalue4 = capfloor.value(valuation_date, libor_curve, model4)
        fvalue5 = capfloor.value(valuation_date, libor_curve, model5)
        fvalue6 = capfloor.value(valuation_date, libor_curve, model6)
        testCases.print("FLR", k, fvalue1, fvalue2, fvalue3, fvalue4, fvalue5, fvalue6)

###############################################################################
# PUT CALL CHECK
###############################################################################

    testCases.header("LABEL", "STRIKE", "BLK", "BLK_SHFTD", "SABR",
                     "SABR SHFTD", "HW", "BACH")

    for k in strikes:
        capFloorType = FinCapFloorTypes.CAP
        capfloor = IborCapFloor(start_date, maturity_date, capFloorType, k)
        cvalue1 = capfloor.value(valuation_date, libor_curve, model1)
        cvalue2 = capfloor.value(valuation_date, libor_curve, model2)
        cvalue3 = capfloor.value(valuation_date, libor_curve, model3)
        cvalue4 = capfloor.value(valuation_date, libor_curve, model4)
        cvalue5 = capfloor.value(valuation_date, libor_curve, model5)
        cvalue6 = capfloor.value(valuation_date, libor_curve, model6)

        capFloorType = FinCapFloorTypes.FLOOR
        capfloor = IborCapFloor(start_date, maturity_date, capFloorType, k)
        fvalue1 = capfloor.value(valuation_date, libor_curve, model1)
        fvalue2 = capfloor.value(valuation_date, libor_curve, model2)
        fvalue3 = capfloor.value(valuation_date, libor_curve, model3)
        fvalue4 = capfloor.value(valuation_date, libor_curve, model4)
        fvalue5 = capfloor.value(valuation_date, libor_curve, model5)
        fvalue6 = capfloor.value(valuation_date, libor_curve, model6)

        pcvalue1 = cvalue1 - fvalue1
        pcvalue2 = cvalue2 - fvalue2
        pcvalue3 = cvalue3 - fvalue3
        pcvalue4 = cvalue4 - fvalue4
        pcvalue5 = cvalue5 - fvalue5
        pcvalue6 = cvalue6 - fvalue6

        testCases.print("PUT_CALL", k, pcvalue1, pcvalue2, pcvalue3,
                        pcvalue4, pcvalue5, pcvalue6)

###############################################################################


def test_FinIborCapFloorVolCurve():
    """ Aim here is to price cap and caplets using cap and caplet vols and to
    demonstrate they are the same - NOT SURE THAT HULLS BOOKS FORMULA WORKS FOR
    OPTIONS. """

    todayDate = Date(20, 6, 2019)
    valuation_date = todayDate
    maturity_date = valuation_date.addTenor("3Y")
    day_count_type = DayCountTypes.THIRTY_E_360
    frequency = FrequencyTypes.ANNUAL

    k = 0.04
    capFloorType = FinCapFloorTypes.CAP
    capFloor = IborCapFloor(valuation_date,
                            maturity_date,
                            capFloorType,
                            k,
                            None,
                            frequency,
                            day_count_type)

    capVolDates = Schedule(valuation_date,
                           valuation_date.addTenor("10Y"),
                           frequency)._generate()

    flat_rate = 0.04
    libor_curve = DiscountCurveFlat(valuation_date,
                                    flat_rate,
                                    frequency,
                                    day_count_type)

    flat = False
    if flat is True:
        capVolatilities = [20.0] * 11
        capVolatilities[0] = 0.0
    else:
        capVolatilities = [0.00, 15.50, 18.25, 17.91, 17.74, 17.27,
                           16.79, 16.30, 16.01, 15.76, 15.54]

    capVolatilities = np.array(capVolatilities)/100.0
    capVolatilities[0] = 0.0

    volCurve = IborCapVolCurve(valuation_date,
                               capVolDates,
                               capVolatilities,
                               day_count_type)

#    print(volCurve._capletGammas)

    # Value cap using a single flat cap volatility
    tcap = (maturity_date - valuation_date) / gDaysInYear
    vol = volCurve.capVol(maturity_date)
    model = FinModelBlack(vol)
    valueCap = capFloor.value(valuation_date, libor_curve, model)
#    print("CAP T", tcap, "VOL:", vol, "VALUE OF CAP:", valueCap)

    # Value cap by breaking it down into caplets using caplet vols
    vCaplets = 0.0
    capletStartDate = capFloor._capFloorLetDates[1]
    testCases.header("START", "END", "VOL", "VALUE")

    for capletEndDate in capFloor._capFloorLetDates[2:]:
        vol = volCurve.capletVol(capletEndDate)
        modelCaplet = FinModelBlack(vol)
        vCaplet = capFloor.valueCapletFloorLet(valuation_date,
                                               capletStartDate,
                                               capletEndDate,
                                               libor_curve,
                                               modelCaplet)

        vCaplets += vCaplet
        testCases.print("%12s" % capletStartDate,
                        "%s" % capletEndDate,
                        "%9.5f" % (vol*100.0),
                        "%9.5f" % vCaplet)

        capletStartDate = capletEndDate

    testCases.header("LABEL", "VALUE")
    testCases.print("CAPLETS->CAP: ", vCaplets)

###############################################################################


def test_FinIborCapletHull():

    #  Hull Page 703, example 29.3
    todayDate = Date(20, 6, 2019)
    valuation_date = todayDate
    maturity_date = valuation_date.addTenor("2Y")
    libor_curve = DiscountCurveFlat(valuation_date,
                                    0.070,
                                    FrequencyTypes.QUARTERLY,
                                    DayCountTypes.THIRTY_E_360)

    k = 0.08
    capFloorType = FinCapFloorTypes.CAP
    capFloor = IborCapFloor(valuation_date,
                            maturity_date,
                            capFloorType,
                            k,
                            None,
                            FrequencyTypes.QUARTERLY,
                            DayCountTypes.THIRTY_E_360)

    # Value cap using a single flat cap volatility
    model = FinModelBlack(0.20)
    capFloor.value(valuation_date, libor_curve, model)

    # Value cap by breaking it down into caplets using caplet vols
    capletStartDate = valuation_date.addTenor("1Y")
    capletEndDate = capletStartDate.addTenor("3M")

    vCaplet = capFloor.valueCapletFloorLet(valuation_date,
                                           capletStartDate,
                                           capletEndDate,
                                           libor_curve,
                                           model)

    # Cannot match Hull due to dates being adjusted
    testCases.header("CORRECT PRICE", "MODEL_PRICE")
    testCases.print(517.29, vCaplet)

###############################################################################


def test_FinIborCapFloorQLExample():

    valuation_date = Date(14, 6, 2016)

    dates = [Date(14, 6, 2016), Date(14, 9, 2016),
             Date(14, 12, 2016), Date(14, 6, 2017),
             Date(14, 6, 2019), Date(14, 6, 2021),
             Date(15, 6, 2026), Date(16, 6, 2031),
             Date(16, 6, 2036), Date(14, 6, 2046)]

    rates = [0.000000, 0.006616, 0.007049, 0.007795,
             0.009599, 0.011203, 0.015068, 0.017583,
             0.018998, 0.020080]

    freq_type = FrequencyTypes.ANNUAL
    day_count_type = DayCountTypes.ACT_ACT_ISDA

    discount_curve = DiscountCurveZeros(valuation_date,
                                        dates,
                                        rates,
                                        freq_type,
                                        day_count_type,
                                        InterpTypes.LINEAR_ZERO_RATES)

    start_date = Date(14, 6, 2016)
    end_date = Date(14, 6, 2026)
    calendar_type = CalendarTypes.UNITED_STATES
    bus_day_adjust_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    freq_type = FrequencyTypes.QUARTERLY
    date_gen_rule_type = DateGenRuleTypes.FORWARD
    lastFixing = 0.0065560
    notional = 1000000
    day_count_type = DayCountTypes.ACT_360
    option_type = FinCapFloorTypes.CAP
    strikeRate = 0.02

    cap = IborCapFloor(start_date, end_date, option_type, strikeRate,
                       lastFixing, freq_type, day_count_type, notional,
                       calendar_type, bus_day_adjust_type, date_gen_rule_type)

    blackVol = 0.547295
    model = FinModelBlack(blackVol)

    start = time.time()
    numRepeats = 10
    for i in range(0, numRepeats):
        v = cap.value(valuation_date, discount_curve, model)

    end = time.time()
    period = end - start
#    print(v, period/numRepeats)

###############################################################################


test_FinIborCapletHull()
test_FinIborCapFloorVolCurve()
test_FinIborCapFloor()
test_FinIborCapFloorQLExample()
testCases.compareTestCases()
