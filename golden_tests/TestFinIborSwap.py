###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys

sys.path.append("..")

import numpy as np
from financepy.utils.math import ONE_MILLION
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.rates.ibor_fra import IborFRA
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
from financepy.utils.global_types import SwapTypes
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.interpolator import InterpTypes
from FinTestCases import FinTestCases, globalTestCaseMode


test_cases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def buildIborSingleCurve(value_dt):

    settle_dt = value_dt.add_days(2)
    dc_type = DayCountTypes.ACT_360

    depos = []
    fras = []
    swaps = []

    maturity_dt = settle_dt.add_months(1)
    depo1 = IborDeposit(value_dt, maturity_dt, -0.00251, dc_type)
    depos.append(depo1)

    # Series of 1M futures
    start_dt = settle_dt.next_imm_date()
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.0023, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00234, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00225, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00226, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00219, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00213, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00186, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00189, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00175, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00143, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00126, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00126, dc_type)
    fras.append(fra)

    ###########################################################################
    ###########################################################################
    ###########################################################################
    ###########################################################################

    fixed_freq = FrequencyTypes.ANNUAL
    dc_type = DayCountTypes.THIRTY_E_360
    fixed_leg_type = SwapTypes.PAY

    #######################################
    maturity_dt = settle_dt.add_months(24)
    swap_rate = -0.001506
    swap1 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap1)

    #######################################
    maturity_dt = settle_dt.add_months(36)
    swap_rate = -0.000185
    swap2 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap2)

    #######################################
    maturity_dt = settle_dt.add_months(48)
    swap_rate = 0.001358
    swap3 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap3)

    #######################################
    maturity_dt = settle_dt.add_months(60)
    swap_rate = 0.0027652
    swap4 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap4)

    #######################################
    maturity_dt = settle_dt.add_months(72)
    swap_rate = 0.0041539
    swap5 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap5)

    #######################################
    maturity_dt = settle_dt.add_months(84)
    swap_rate = 0.0054604
    swap6 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap6)

    #######################################
    maturity_dt = settle_dt.add_months(96)
    swap_rate = 0.006674
    swap7 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap7)

    #######################################
    maturity_dt = settle_dt.add_months(108)
    swap_rate = 0.007826
    swap8 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap8)

    #######################################
    maturity_dt = settle_dt.add_months(120)
    swap_rate = 0.008821
    swap9 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap9)

    #######################################
    maturity_dt = settle_dt.add_months(132)
    swap_rate = 0.0097379
    swap10 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap10)

    #######################################
    maturity_dt = settle_dt.add_months(144)
    swap_rate = 0.0105406
    swap11 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap11)

    #######################################
    maturity_dt = settle_dt.add_months(180)
    swap_rate = 0.0123927
    swap12 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap12)

    #######################################
    maturity_dt = settle_dt.add_months(240)
    swap_rate = 0.0139882
    swap13 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap13)

    #######################################
    maturity_dt = settle_dt.add_months(300)
    swap_rate = 0.0144972
    swap14 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap14)

    #######################################
    maturity_dt = settle_dt.add_months(360)
    swap_rate = 0.0146081
    swap15 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap15)

    #######################################
    maturity_dt = settle_dt.add_months(420)
    swap_rate = 0.01461897
    swap16 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap16)

    #######################################
    maturity_dt = settle_dt.add_months(480)
    swap_rate = 0.014567455
    swap17 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap17)

    #######################################
    maturity_dt = settle_dt.add_months(540)
    swap_rate = 0.0140826
    swap18 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap18)

    #######################################
    maturity_dt = settle_dt.add_months(600)
    swap_rate = 0.01436822
    swap19 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap19)

    ########################################

    libor_curve = IborSingleCurve(value_dt, depos, fras, swaps)

    test_cases.header("LABEL", "DATE", "VALUE")

    """ Check calibration """
    for depo in depos:
        v = depo.value(settle_dt, libor_curve)
        test_cases.print("DEPO VALUE:", depo.maturity_dt, v)

    for fra in fras:
        v = fra.value(settle_dt, libor_curve)
        test_cases.print("FRA VALUE:", fra.maturity_dt, v)

    for swap in swaps:
        v = swap.value(settle_dt, libor_curve)
        test_cases.print("SWAP VALUE:", swap.maturity_dt, v)

    return libor_curve


###############################################################################


def test_LiborSwap():

    # I have tried to reproduce the example from the blog by Ioannis Rigopoulos
    # https://blog.deriscope.com/index.php/en/excel-interest-rate-swap-price-dual-bootstrapping-curve
    start_dt = Date(27, 12, 2017)
    end_dt = Date(27, 12, 2067)

    fixed_cpn = 0.015
    fixed_freq_type = FrequencyTypes.ANNUAL
    fixed_dc_type = DayCountTypes.THIRTY_E_360

    float_spread = 0.0
    float_freq_type = FrequencyTypes.SEMI_ANNUAL
    float_dc_type = DayCountTypes.ACT_360
    firstFixing = -0.00268

    swap_cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    fixed_leg_type = SwapTypes.RECEIVE

    notional = 10.0 * ONE_MILLION

    swap = IborSwap(
        start_dt,
        end_dt,
        fixed_leg_type,
        fixed_cpn,
        fixed_freq_type,
        fixed_dc_type,
        notional,
        float_spread,
        float_freq_type,
        float_dc_type,
        swap_cal_type,
        bd_type,
        dg_type,
    )

    """ Now perform a valuation after the swap has seasoned but with the
    same curve being used for discounting and working out the implied
    future Libor rates. """

    value_dt = Date(30, 11, 2018)
    settle_dt = value_dt.add_days(2)
    libor_curve = buildIborSingleCurve(value_dt)
    v = swap.value(settle_dt, libor_curve, libor_curve, firstFixing)

    v_bbg = 388147.0
    test_cases.header("LABEL", "VALUE")
    test_cases.print("SWAP_VALUE USING ONE_CURVE", v)
    test_cases.print("BLOOMBERG VALUE", v_bbg)
    test_cases.print("DIFFERENCE VALUE", v_bbg - v)


###############################################################################


def test_dp_example():

    #  http://www.derivativepricing.com/blogpage.asp?id=8

    start_dt = Date(14, 11, 2011)
    end_dt = Date(14, 11, 2016)
    fixed_freq_type = FrequencyTypes.SEMI_ANNUAL
    swap_cal_type = CalendarTypes.TARGET
    bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    fixed_dc_type = DayCountTypes.THIRTY_E_360_ISDA
    fixed_leg_type = SwapTypes.PAY
    fixed_cpn = 0.0124
    notional = ONE_MILLION

    swap = IborSwap(
        start_dt,
        end_dt,
        fixed_leg_type,
        fixed_cpn=fixed_cpn,
        fixed_freq_type=fixed_freq_type,
        fixed_dc_type=fixed_dc_type,
        float_freq_type=FrequencyTypes.SEMI_ANNUAL,
        float_dc_type=DayCountTypes.ACT_360,
        notional=notional,
        cal_type=swap_cal_type,
        bd_type=bd_type,
        dg_type=dg_type,
    )

    #    swap.printFixedLegFlows()

    dts = [
        Date(14, 11, 2011),
        Date(14, 5, 2012),
        Date(14, 11, 2012),
        Date(14, 5, 2013),
        Date(14, 11, 2013),
        Date(14, 5, 2014),
        Date(14, 11, 2014),
        Date(14, 5, 2015),
        Date(16, 11, 2015),
        Date(16, 5, 2016),
        Date(14, 11, 2016),
    ]

    dfs = [
        0.9999843,
        0.9966889,
        0.9942107,
        0.9911884,
        0.9880738,
        0.9836490,
        0.9786276,
        0.9710461,
        0.9621778,
        0.9514315,
        0.9394919,
    ]

    value_dt = start_dt

    curve = DiscountCurve(
        value_dt, dts, np.array(dfs), InterpTypes.FLAT_FWD_RATES
    )

    v = swap.value(value_dt, curve, curve)

    #    swap.print_fixed_leg_pv()
    #    swap.print_float_leg_pv()

    # This is essentially zero
    test_cases.header("LABEL", "VALUE")
    test_cases.print("Swap Value on a Notional of $1M:", v)


###############################################################################


test_LiborSwap()
test_dp_example()
test_cases.compareTestCases()
