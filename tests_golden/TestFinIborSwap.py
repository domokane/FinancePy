###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

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
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def buildIborSingleCurve(valuation_date):

    settlement_date = valuation_date.add_days(2)
    dcType = DayCountTypes.ACT_360

    depos = []
    fras = []
    swaps = []

    maturity_date = settlement_date.add_months(1)
    depo1 = IborDeposit(valuation_date, maturity_date, -0.00251, dcType)
    depos.append(depo1)

    # Series of 1M futures
    start_date = settlement_date.next_imm_date()
    end_date = start_date.add_months(1)
    fra = IborFRA(start_date, end_date, -0.0023, dcType)
    fras.append(fra)

    start_date = start_date.add_months(1)
    end_date = start_date.add_months(1)
    fra = IborFRA(start_date, end_date, -0.00234, dcType)
    fras.append(fra)

    start_date = start_date.add_months(1)
    end_date = start_date.add_months(1)
    fra = IborFRA(start_date, end_date, -0.00225, dcType)
    fras.append(fra)

    start_date = start_date.add_months(1)
    end_date = start_date.add_months(1)
    fra = IborFRA(start_date, end_date, -0.00226, dcType)
    fras.append(fra)

    start_date = start_date.add_months(1)
    end_date = start_date.add_months(1)
    fra = IborFRA(start_date, end_date, -0.00219, dcType)
    fras.append(fra)

    start_date = start_date.add_months(1)
    end_date = start_date.add_months(1)
    fra = IborFRA(start_date, end_date, -0.00213, dcType)
    fras.append(fra)

    start_date = start_date.add_months(1)
    end_date = start_date.add_months(1)
    fra = IborFRA(start_date, end_date, -0.00186, dcType)
    fras.append(fra)

    start_date = start_date.add_months(1)
    end_date = start_date.add_months(1)
    fra = IborFRA(start_date, end_date, -0.00189, dcType)
    fras.append(fra)

    start_date = start_date.add_months(1)
    end_date = start_date.add_months(1)
    fra = IborFRA(start_date, end_date, -0.00175, dcType)
    fras.append(fra)

    start_date = start_date.add_months(1)
    end_date = start_date.add_months(1)
    fra = IborFRA(start_date, end_date, -0.00143, dcType)
    fras.append(fra)

    start_date = start_date.add_months(1)
    end_date = start_date.add_months(1)
    fra = IborFRA(start_date, end_date, -0.00126, dcType)
    fras.append(fra)

    start_date = start_date.add_months(1)
    end_date = start_date.add_months(1)
    fra = IborFRA(start_date, end_date, -0.00126, dcType)
    fras.append(fra)

    ###########################################################################
    ###########################################################################
    ###########################################################################
    ###########################################################################

    fixedFreq = FrequencyTypes.ANNUAL
    dcType = DayCountTypes.THIRTY_E_360
    fixed_leg_type = SwapTypes.PAY

    #######################################
    maturity_date = settlement_date.add_months(24)
    swap_rate = -0.001506
    swap1 = IborSwap(settlement_date, maturity_date, fixed_leg_type,
                     swap_rate, fixedFreq, dcType)
    swaps.append(swap1)

    #######################################
    maturity_date = settlement_date.add_months(36)
    swap_rate = -0.000185
    swap2 = IborSwap(settlement_date, maturity_date, fixed_leg_type,
                     swap_rate, fixedFreq, dcType)
    swaps.append(swap2)

    #######################################
    maturity_date = settlement_date.add_months(48)
    swap_rate = 0.001358
    swap3 = IborSwap(settlement_date, maturity_date, fixed_leg_type,
                     swap_rate, fixedFreq, dcType)
    swaps.append(swap3)

    #######################################
    maturity_date = settlement_date.add_months(60)
    swap_rate = 0.0027652
    swap4 = IborSwap(settlement_date, maturity_date, fixed_leg_type,
                     swap_rate, fixedFreq, dcType)
    swaps.append(swap4)

    #######################################
    maturity_date = settlement_date.add_months(72)
    swap_rate = 0.0041539
    swap5 = IborSwap(settlement_date, maturity_date, fixed_leg_type,
                     swap_rate, fixedFreq, dcType)
    swaps.append(swap5)

    #######################################
    maturity_date = settlement_date.add_months(84)
    swap_rate = 0.0054604
    swap6 = IborSwap(settlement_date, maturity_date, fixed_leg_type,
                     swap_rate, fixedFreq, dcType)
    swaps.append(swap6)

    #######################################
    maturity_date = settlement_date.add_months(96)
    swap_rate = 0.006674
    swap7 = IborSwap(settlement_date, maturity_date, fixed_leg_type,
                     swap_rate, fixedFreq, dcType)
    swaps.append(swap7)

    #######################################
    maturity_date = settlement_date.add_months(108)
    swap_rate = 0.007826
    swap8 = IborSwap(settlement_date, maturity_date, fixed_leg_type,
                     swap_rate, fixedFreq, dcType)
    swaps.append(swap8)

    #######################################
    maturity_date = settlement_date.add_months(120)
    swap_rate = 0.008821
    swap9 = IborSwap(settlement_date, maturity_date, fixed_leg_type,
                     swap_rate, fixedFreq, dcType)
    swaps.append(swap9)

    #######################################
    maturity_date = settlement_date.add_months(132)
    swap_rate = 0.0097379
    swap10 = IborSwap(settlement_date, maturity_date, fixed_leg_type,
                      swap_rate, fixedFreq, dcType)
    swaps.append(swap10)

    #######################################
    maturity_date = settlement_date.add_months(144)
    swap_rate = 0.0105406
    swap11 = IborSwap(settlement_date, maturity_date, fixed_leg_type,
                      swap_rate, fixedFreq, dcType)
    swaps.append(swap11)

    #######################################
    maturity_date = settlement_date.add_months(180)
    swap_rate = 0.0123927
    swap12 = IborSwap(settlement_date, maturity_date, fixed_leg_type,
                      swap_rate, fixedFreq, dcType)
    swaps.append(swap12)

    #######################################
    maturity_date = settlement_date.add_months(240)
    swap_rate = 0.0139882
    swap13 = IborSwap(settlement_date, maturity_date, fixed_leg_type,
                      swap_rate, fixedFreq, dcType)
    swaps.append(swap13)

    #######################################
    maturity_date = settlement_date.add_months(300)
    swap_rate = 0.0144972
    swap14 = IborSwap(settlement_date, maturity_date, fixed_leg_type,
                      swap_rate, fixedFreq, dcType)
    swaps.append(swap14)

    #######################################
    maturity_date = settlement_date.add_months(360)
    swap_rate = 0.0146081
    swap15 = IborSwap(settlement_date, maturity_date, fixed_leg_type,
                      swap_rate, fixedFreq, dcType)
    swaps.append(swap15)

    #######################################
    maturity_date = settlement_date.add_months(420)
    swap_rate = 0.01461897
    swap16 = IborSwap(settlement_date, maturity_date, fixed_leg_type,
                      swap_rate, fixedFreq, dcType)
    swaps.append(swap16)

    #######################################
    maturity_date = settlement_date.add_months(480)
    swap_rate = 0.014567455
    swap17 = IborSwap(settlement_date, maturity_date, fixed_leg_type,
                      swap_rate, fixedFreq, dcType)
    swaps.append(swap17)

    #######################################
    maturity_date = settlement_date.add_months(540)
    swap_rate = 0.0140826
    swap18 = IborSwap(settlement_date, maturity_date, fixed_leg_type,
                      swap_rate, fixedFreq, dcType)
    swaps.append(swap18)

    #######################################
    maturity_date = settlement_date.add_months(600)
    swap_rate = 0.01436822
    swap19 = IborSwap(settlement_date, maturity_date, fixed_leg_type,
                      swap_rate, fixedFreq, dcType)
    swaps.append(swap19)

    ########################################

    libor_curve = IborSingleCurve(valuation_date, depos, fras, swaps)

    testCases.header("LABEL", "DATE", "VALUE")

    """ Check calibration """
    for depo in depos:
        v = depo.value(settlement_date, libor_curve)
        testCases.print("DEPO VALUE:", depo._maturity_date, v)

    for fra in fras:
        v = fra.value(settlement_date, libor_curve)
        testCases.print("FRA VALUE:", fra._maturity_date, v)

    for swap in swaps:
        v = swap.value(settlement_date, libor_curve)
        testCases.print("SWAP VALUE:", swap._maturity_date, v)

    return libor_curve

###############################################################################


def test_LiborSwap():

    # I have tried to reproduce the example from the blog by Ioannis Rigopoulos
    # https://blog.deriscope.com/index.php/en/excel-interest-rate-swap-price-dual-bootstrapping-curve
    start_date = Date(27, 12, 2017)
    end_date = Date(27, 12, 2067)

    fixed_coupon = 0.015
    fixedFreqType = FrequencyTypes.ANNUAL
    fixed_day_count_type = DayCountTypes.THIRTY_E_360

    float_spread = 0.0
    floatFreqType = FrequencyTypes.SEMI_ANNUAL
    float_day_count_type = DayCountTypes.ACT_360
    firstFixing = -0.00268

    swapCalendarType = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
    fixed_leg_type = SwapTypes.RECEIVE

    notional = 10.0 * ONE_MILLION

    swap = IborSwap(start_date,
                    end_date,
                    fixed_leg_type,
                    fixed_coupon,
                    fixedFreqType,
                    fixed_day_count_type,
                    notional,
                    float_spread,
                    floatFreqType,
                    float_day_count_type,
                    swapCalendarType,
                    bus_day_adjust_type,
                    date_gen_rule_type)

    """ Now perform a valuation after the swap has seasoned but with the
    same curve being used for discounting and working out the implied
    future Libor rates. """

    valuation_date = Date(30, 11, 2018)
    settlement_date = valuation_date.add_days(2)
    libor_curve = buildIborSingleCurve(valuation_date)
    v = swap.value(settlement_date, libor_curve, libor_curve, firstFixing)

    v_bbg = 388147.0
    testCases.header("LABEL", "VALUE")
    testCases.print("SWAP_VALUE USING ONE_CURVE", v)
    testCases.print("BLOOMBERG VALUE", v_bbg)
    testCases.print("DIFFERENCE VALUE", v_bbg - v)

###############################################################################


def test_dp_example():

    #  http://www.derivativepricing.com/blogpage.asp?id=8

    start_date = Date(14, 11, 2011)
    end_date = Date(14, 11, 2016)
    fixedFreqType = FrequencyTypes.SEMI_ANNUAL
    swapCalendarType = CalendarTypes.TARGET
    bus_day_adjust_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
    fixed_day_count_type = DayCountTypes.THIRTY_E_360_ISDA
    fixed_leg_type = SwapTypes.PAY
    fixed_coupon = 0.0124
    notional = ONE_MILLION

    swap = IborSwap(start_date,
                    end_date,
                    fixed_leg_type,
                    fixed_coupon=fixed_coupon,
                    fixed_freq_type=fixedFreqType,
                    fixed_day_count_type=fixed_day_count_type,
                    float_freq_type=FrequencyTypes.SEMI_ANNUAL,
                    float_day_count_type=DayCountTypes.ACT_360,
                    notional=notional,
                    calendar_type=swapCalendarType,
                    bus_day_adjust_type=bus_day_adjust_type,
                    date_gen_rule_type=date_gen_rule_type)

#    swap.printFixedLegFlows()

    dts = [Date(14, 11, 2011), Date(14, 5, 2012), Date(14, 11, 2012),
           Date(14, 5, 2013), Date(14, 11, 2013), Date(14, 5, 2014),
           Date(14, 11, 2014), Date(14, 5, 2015), Date(16, 11, 2015),
           Date(16, 5, 2016), Date(14, 11, 2016)]

    dfs = [0.9999843, 0.9966889, 0.9942107, 0.9911884, 0.9880738, 0.9836490,
           0.9786276, 0.9710461, 0.9621778, 0.9514315, 0.9394919]

    valuation_date = start_date

    curve = DiscountCurve(valuation_date, dts, np.array(dfs),
                          InterpTypes.FLAT_FWD_RATES)

    v = swap.value(valuation_date, curve, curve)

#    swap.print_fixed_leg_pv()
#    swap.print_float_leg_pv()

    # This is essentially zero
    testCases.header("LABEL", "VALUE")
    testCases.print("Swap Value on a Notional of $1M:", v)

###############################################################################


test_LiborSwap()
test_dp_example()
testCases.compareTestCases()
