###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import matplotlib.pyplot as plt
import numpy as np
import time as time

import sys
sys.path.append("..")

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.calendar import CalendarTypes
from financepy.products.rates.FinIborSingleCurveOLD import FinIborSingleCurveOLD
from financepy.products.rates.ibor_fra import IborFRA
from financepy.products.rates.ibor_future import IborFuture
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.FinIborSwapOLD import FinIborSwapOLD
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.market.curves.interpolator import InterpTypes
from financepy.utils.math import ONE_MILLION
from financepy.utils.global_types import FinSwapTypes

testCases = FinTestCases(__file__, globalTestCaseMode)

PLOT_GRAPHS = False

###############################################################################


def test_FinIborDepositsOnly():

    # I have used the following useful blog post by Ioannis Rigopoulos for this
    # https://blog.deriscope.com/index.php/en/yield-curve-excel-quantlib-deposit

    valuation_date = Date(23, 2, 2018)

    spot_days = 0
    settlement_date = valuation_date.addWeekDays(spot_days)

    depoDCCType = DayCountTypes.ACT_360
    notional = 100.0
    calendar_type = CalendarTypes.TARGET
    depos = []

    # 1 month
    deposit_rate = 0.04
    maturity_date = settlement_date.addMonths(1)
    depo = IborDeposit(settlement_date, maturity_date, deposit_rate,
                       depoDCCType, notional, calendar_type)
    depos.append(depo)

    # 2 months
    deposit_rate = 0.04
    maturity_date = settlement_date.addMonths(2)
    depo = IborDeposit(settlement_date, maturity_date, deposit_rate,
                       depoDCCType, notional, calendar_type)
    depos.append(depo)

    # 6 months
    deposit_rate = 0.04
    maturity_date = settlement_date.addMonths(6)
    depo = IborDeposit(settlement_date, maturity_date, deposit_rate,
                       depoDCCType, notional, calendar_type)
    depos.append(depo)

    # 1 year
    deposit_rate = 0.04
    maturity_date = settlement_date.addMonths(12)
    depo = IborDeposit(settlement_date, maturity_date, deposit_rate,
                       depoDCCType, notional, calendar_type)
    depos.append(depo)

    fras = []
    swaps = []

    libor_curve = FinIborSingleCurveOLD(valuation_date,
                                       depos,
                                       fras,
                                       swaps)

    testCases.header("LABEL", "DATE", "VALUE")

    """ Check calibration """
    for depo in depos:
        v = depo.value(settlement_date, libor_curve)
        testCases.print("DEPO", depo._maturity_date, v)

###############################################################################


def test_FinIborFRAsOnly():

    # TO DO FIX THIS
    valuation_date = Date(23, 2, 2018)

    spot_days = 0
    settlement_date = valuation_date.addWeekDays(spot_days)

    depoDCCType = DayCountTypes.ACT_360
    notional = 100.0

    payFixed = True

    calendar_type = CalendarTypes.TARGET
    fras = []

    # 1 x 4 FRA
    fraRate = 0.04
    fraSettlementDate = settlement_date.addMonths(1)
    fraMaturityDate = settlement_date.addMonths(4)
    fra = IborFRA(fraSettlementDate, fraMaturityDate, fraRate,
                  depoDCCType, notional, payFixed, calendar_type)
    fras.append(fra)

    # 4 x 7 FRA
    fraRate = 0.08
    fraSettlementDate = settlement_date.addMonths(4)
    fraMaturityDate = settlement_date.addMonths(7)
    fra = IborFRA(fraSettlementDate, fraMaturityDate, fraRate,
                  depoDCCType, notional, payFixed, calendar_type)
    fras.append(fra)

    depos = []
    swaps = []

    libor_curve = FinIborSingleCurveOLD(valuation_date,
                                       depos,
                                       fras,
                                       swaps)

    testCases.header("DATE", "MATDATE", "VALUE")

    """ Check calibration """
    for fra in fras:
        v = fra.value(settlement_date, libor_curve)
        testCases.print("FRA:", fra._maturity_date, v)

###############################################################################


def test_FinIborDepositsFRAsSwaps():

    valuation_date = Date(18, 9, 2019)

    dccType = DayCountTypes.THIRTY_E_360_ISDA
    depos = []

    spot_days = 0
    settlement_date = valuation_date.addWeekDays(spot_days)

    deposit_rate = 0.050
    maturity_date = settlement_date.addMonths(1)
    depo = IborDeposit(settlement_date, maturity_date, deposit_rate, dccType)
    depos.append(depo)

    maturity_date = settlement_date.addMonths(2)
    depo = IborDeposit(settlement_date, maturity_date, deposit_rate, dccType)
    depos.append(depo)

    maturity_date = settlement_date.addMonths(3)
    depo = IborDeposit(settlement_date, maturity_date, deposit_rate, dccType)
    depos.append(depo)

    maturity_date = settlement_date.addMonths(6)
    depo = IborDeposit(settlement_date, maturity_date, deposit_rate, dccType)
    depos.append(depo)

    maturity_date = settlement_date.addMonths(9)
    depo = IborDeposit(settlement_date, maturity_date, deposit_rate, dccType)
    depos.append(depo)

    maturity_date = settlement_date.addMonths(12)
    depo = IborDeposit(settlement_date, maturity_date, deposit_rate, dccType)
    depos.append(depo)

    fras = []
    # 1 x 4 FRA
    fraRate = 0.04
    fraSettlementDate = settlement_date.addMonths(9)
    fraMaturityDate = settlement_date.addMonths(13)
    fra = IborFRA(fraSettlementDate, fraMaturityDate, fraRate, dccType)
    fras.append(fra)

    # 4 x 7 FRA
    fraRate = 0.03
    fraSettlementDate = settlement_date.addMonths(13)
    fraMaturityDate = settlement_date.addMonths(17)
    fra = IborFRA(fraSettlementDate, fraMaturityDate, fraRate, dccType)
    fras.append(fra)

    # 4 x 7 FRA
    fraRate = 0.07
    fraSettlementDate = settlement_date.addMonths(17)
    fraMaturityDate = settlement_date.addMonths(21)
    fra = IborFRA(fraSettlementDate, fraMaturityDate, fraRate, dccType)
    fras.append(fra)

    swaps = []
    fixedDCCType = DayCountTypes.ACT_365F
    fixedFreqType = FrequencyTypes.SEMI_ANNUAL

    swap_rate = 0.05
#    maturity_date = settlement_date.addMonths(24)
#    swap = FinIborSwapOLD(settlement_date, maturity_date, swap_rate, fixedFreqType,
#                        fixedDCCType)
#    swaps.append(swap)

    fixed_leg_type = FinSwapTypes.PAY
    maturity_date = settlement_date.addMonths(36)
    swap = FinIborSwapOLD(settlement_date, maturity_date, fixed_leg_type, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settlement_date.addMonths(48)
    swap = FinIborSwapOLD(settlement_date, maturity_date, fixed_leg_type, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settlement_date.addMonths(60)
    swap = FinIborSwapOLD(settlement_date, maturity_date, fixed_leg_type, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settlement_date.addMonths(72)
    swap = FinIborSwapOLD(settlement_date, maturity_date, fixed_leg_type, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settlement_date.addMonths(84)
    swap = FinIborSwapOLD(settlement_date, maturity_date, fixed_leg_type, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settlement_date.addMonths(96)
    swap = FinIborSwapOLD(settlement_date, maturity_date, fixed_leg_type, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settlement_date.addMonths(108)
    swap = FinIborSwapOLD(settlement_date, maturity_date, fixed_leg_type, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settlement_date.addMonths(120)
    swap = FinIborSwapOLD(settlement_date, maturity_date, fixed_leg_type, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settlement_date.addMonths(132)
    swap = FinIborSwapOLD(settlement_date, maturity_date, fixed_leg_type, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settlement_date.addMonths(144)
    swap = FinIborSwapOLD(settlement_date, maturity_date, fixed_leg_type, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settlement_date.addMonths(180)
    swap = FinIborSwapOLD(settlement_date, maturity_date, fixed_leg_type, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settlement_date.addMonths(240)
    swap = FinIborSwapOLD(settlement_date, maturity_date, fixed_leg_type, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settlement_date.addMonths(300)
    swap = FinIborSwapOLD(settlement_date, maturity_date, fixed_leg_type, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    maturity_date = settlement_date.addMonths(360)
    swap = FinIborSwapOLD(settlement_date, maturity_date, fixed_leg_type, swap_rate,
                        fixedFreqType,
                        fixedDCCType)
    swaps.append(swap)

    libor_curve = FinIborSingleCurveOLD(valuation_date,
                                       depos,
                                       fras,
                                       swaps)

    df = libor_curve.df(settlement_date)

    testCases.header("SETTLEMENT DATE", "DF")
    testCases.print(str(settlement_date), df)
    testCases.header("DATE", "DF")

    for deposit in depos:
        df = libor_curve.df(deposit._maturity_date)
        testCases.print(str(deposit._maturity_date), df)

    for swap in swaps:
        df = libor_curve.df(swap._maturity_date)
        testCases.print(str(swap._maturity_date), df)


###############################################################################
###############################################################################
# https://github.com/vilen22/curve-building/blob/master/Bloomberg%20Curve%20Building%20Replication.xlsx
###############################################################################
###############################################################################
# AGREEEMENT IS VERY CLOSE - NOT SURE ABOUT SIZE OF LAST PAYMENT ON FIXED LEG!


def futureToFRARate(price, convexity):
    futRate = (100-price)/100
    if convexity < 0:
        fraRate = futRate + convexity/100.0
    else:
        fraRate = futRate - convexity/100.0

    return fraRate

###############################################################################


def test_FinIborDepositsFuturesSwaps():

    spotDate = Date(6, 6, 2018)
    spot_days = 0
    settlement_date = spotDate.addWeekDays(spot_days)
    depoDCCType = DayCountTypes.ACT_360
    depos = []
    deposit_rate = 0.0231381
    depo = IborDeposit(settlement_date, "3M", deposit_rate, depoDCCType)
    depos.append(depo)

    deposit_rate = 0.027
    depo = IborDeposit(settlement_date, "3M", deposit_rate, depoDCCType)
    depos.append(depo)

    depos = []
    depo = IborDeposit(settlement_date, "1M", 0.0230, depoDCCType)
    depos.append(depo)
    depo = IborDeposit(settlement_date, "2M", 0.0235, depoDCCType)
    depos.append(depo)
    depo = IborDeposit(settlement_date, "3M", 0.0240, depoDCCType)
    depos.append(depo)

    fras = []

    fraRate = futureToFRARate(97.6675, -0.00005)
    fraSettlementDate = spotDate.nextIMMDate()
    fraMaturityDate = fraSettlementDate.nextIMMDate()
    fra = IborFRA(fraSettlementDate, fraMaturityDate, fraRate, depoDCCType)
    fras.append(fra)

    fraRate = futureToFRARate(97.5200, -0.00060)
    fraSettlementDate = fraMaturityDate
    fraMaturityDate = fraSettlementDate.nextIMMDate()
    fra = IborFRA(fraSettlementDate, fraMaturityDate, fraRate, depoDCCType)
    fras.append(fra)

    fraRate = futureToFRARate(97.3550, -0.00146)
    fraSettlementDate = fraMaturityDate
    fraMaturityDate = fraSettlementDate.nextIMMDate()
    fra = IborFRA(fraSettlementDate, fraMaturityDate, fraRate, depoDCCType)
    fras.append(fra)

    fraRate = futureToFRARate(97.2450, -0.00263)
    fraSettlementDate = fraMaturityDate
    fraMaturityDate = fraSettlementDate.nextIMMDate()
    fra = IborFRA(fraSettlementDate, fraMaturityDate, fraRate, depoDCCType)
    fras.append(fra)

    fraRate = futureToFRARate(97.1450, -0.00411)
    fraSettlementDate = fraMaturityDate
    fraMaturityDate = fraSettlementDate.nextIMMDate()
    fra = IborFRA(fraSettlementDate, fraMaturityDate, fraRate, depoDCCType)
    fras.append(fra)

    fraRate = futureToFRARate(97.0750, -0.00589)
    fraSettlementDate = fraSettlementDate.nextIMMDate()
    fraMaturityDate = fraSettlementDate.nextIMMDate()
    fra = IborFRA(fraSettlementDate, fraMaturityDate, fraRate, depoDCCType)
    fras.append(fra)

    ###########################################################################

    spot_days = 2
    start_date = spotDate.addWeekDays(spot_days)

    swaps = []
    fixed_leg_type = FinSwapTypes.PAY
    fixedDCCType = DayCountTypes.THIRTY_E_360
    fixedFreqType = FrequencyTypes.SEMI_ANNUAL
    floatFreqType = FrequencyTypes.QUARTERLY
    notional = 1000000
    principal = 0.0
    float_spread = 0.0
    floatDCCType = DayCountTypes.ACT_360
    calendar_type = CalendarTypes.US
    busDayAdjustRule = BusDayAdjustTypes.PRECEDING

    swap_rate = 0.02776305

    swap = FinIborSwapOLD(start_date, "2Y", fixed_leg_type, swap_rate,
                        fixedFreqType, fixedDCCType, notional,
                        float_spread, floatFreqType, floatDCCType,
                        calendar_type, busDayAdjustRule)

    swaps.append(swap)

    libor_curve = FinIborSingleCurveOLD(spotDate, depos, fras, swaps)

    times = np.linspace(0.0, 2.0, 25)
    dates = spotDate.addYears(times)
    zero_rates = libor_curve.zero_rate(dates)
    fwd_rates = libor_curve.fwd(dates)

    if PLOT_GRAPHS:
        plt.figure(figsize=(8, 6))
        plt.plot(times, zero_rates*100, label="zero rates")
        plt.plot(times, fwd_rates*100, label="fwd rates")
        plt.xlabel("Times")
        plt.ylabel("CC forward rates")
        plt.legend()

        print("==============================================================")
        for fra in fras:
            print(fra)
        print("==============================================================")

        end_date = spotDate
        df = libor_curve.df(end_date)
        print(end_date, df)

        end_date = settlement_date
        df = libor_curve.df(end_date)
        print(end_date, df)

        end_date = Date(20, 6, 2018)
        df = libor_curve.df(end_date)
        print(end_date, df)

        for depo in depos:
            end_date = depo._maturity_date
            df = libor_curve.df(end_date)
            print(end_date, df)

        for fra in fras:
            end_date = fra._maturity_date
            df = libor_curve.df(end_date)
            print(end_date, df)

        for swap in swaps:
            end_date = swap._maturity_date
            df = libor_curve.df(end_date)
            print(end_date, df)

        swap.print_fixed_leg_pv(spotDate)
        swap.print_float_leg_pv(spotDate)

#        print(libor_curve)

###############################################################################


def test_derivativePricingExample():

    valuation_date = Date(10, 11, 2011)

    dccType = DayCountTypes.ACT_360
    depos = []

    # We do the O/N rate which settles on trade date
    spot_days = 0
    settlement_date = valuation_date.addWeekDays(spot_days)

    deposit_rate = 0.001410
    depo = IborDeposit(settlement_date, "ON", deposit_rate, dccType)
    depos.append(depo)

    spot_days = 1
    settlement_date = valuation_date.addWeekDays(spot_days)

    deposit_rate = 0.001410
    depo = IborDeposit(settlement_date, "TN", deposit_rate, dccType)
    depos.append(depo)

    spot_days = 2
    settlement_date = valuation_date.addWeekDays(spot_days)

    deposit_rate = 0.001910
    depo = IborDeposit(settlement_date, "1W", deposit_rate, dccType)
    depos.append(depo)

    deposit_rate = 0.002090
    depo = IborDeposit(settlement_date, "2W", deposit_rate, dccType)
    depos.append(depo)

    deposit_rate = 0.002490
    depo = IborDeposit(settlement_date, "1M", deposit_rate, dccType)
    depos.append(depo)

    deposit_rate = 0.003450
    depo = IborDeposit(settlement_date, "2M", deposit_rate, dccType)
    depos.append(depo)

    deposit_rate = 0.004570
    depo = IborDeposit(settlement_date, "3M", deposit_rate, dccType)
    depos.append(depo)

    deposit_rate = 0.005230
    depo = IborDeposit(settlement_date, "4M", deposit_rate, dccType)
    depos.append(depo)

    deposit_rate = 0.005860
    depo = IborDeposit(settlement_date, "5M", deposit_rate, dccType)
    depos.append(depo)

    deposit_rate = 0.006540
    depo = IborDeposit(settlement_date, "6M", deposit_rate, dccType)
    depos.append(depo)

    deposit_rate = 0.007080
    depo = IborDeposit(settlement_date, "7M", deposit_rate, dccType)
    depos.append(depo)

    deposit_rate = 0.007540
    depo = IborDeposit(settlement_date, "8M", deposit_rate, dccType)
    depos.append(depo)

    deposit_rate = 0.008080
    depo = IborDeposit(settlement_date, "9M", deposit_rate, dccType)
    depos.append(depo)

    deposit_rate = 0.008570
    depo = IborDeposit(settlement_date, "10M", deposit_rate, dccType)
    depos.append(depo)

    deposit_rate = 0.009130
    depo = IborDeposit(settlement_date, "11M", deposit_rate, dccType)
    depos.append(depo)

    fras = []

    swaps = []
    day_count_type = DayCountTypes.THIRTY_E_360_ISDA
#    day_count_type = DayCountTypes.ACT_360
    freq_type = FrequencyTypes.SEMI_ANNUAL
    fixed_leg_type = FinSwapTypes.PAY
    
    swap_rate = 0.0058
    swap = FinIborSwapOLD(settlement_date, "1Y", fixed_leg_type, swap_rate, freq_type, day_count_type)
    swaps.append(swap)

    swap_rate = 0.0060
    swap = FinIborSwapOLD(settlement_date, "2Y", fixed_leg_type, swap_rate, freq_type, day_count_type)
    swaps.append(swap)

    swap_rate = 0.0072
    swap = FinIborSwapOLD(settlement_date, "3Y", fixed_leg_type, swap_rate, freq_type, day_count_type)
    swaps.append(swap)

    swap_rate = 0.0096
    swap = FinIborSwapOLD(settlement_date, "4Y", fixed_leg_type, swap_rate, freq_type, day_count_type)
    swaps.append(swap)

    swap_rate = 0.0124
    swap = FinIborSwapOLD(settlement_date, "5Y", fixed_leg_type, swap_rate, freq_type, day_count_type)
    swaps.append(swap)

    swap_rate = 0.0173
    swap = FinIborSwapOLD(settlement_date, "7Y", fixed_leg_type, swap_rate, freq_type, day_count_type)
    swaps.append(swap)

    swap_rate = 0.0219
    swap = FinIborSwapOLD(settlement_date, "10Y", fixed_leg_type, swap_rate, freq_type, day_count_type)
    swaps.append(swap)

    swap_rate = 0.0283
    swap = FinIborSwapOLD(settlement_date, "30Y", fixed_leg_type, swap_rate, freq_type, day_count_type)
    swaps.append(swap)

    numRepeats = 10
    start = time.time()

    for _ in range(0, numRepeats):
        _ = FinIborSingleCurveOLD(valuation_date, depos, fras, swaps,
                                  InterpTypes.FLAT_FWD_RATES)

    end = time.time()
    elapsed1 = end - start

    start = time.time()

    for _ in range(0, numRepeats):
        _ = FinIborSingleCurveOLD(valuation_date, depos, fras, swaps,
                                  InterpTypes.FLAT_FWD_RATES)

    end = time.time()
    elapsed2 = end - start

    testCases.header("METHOD", "TIME")
    testCases.print("NON-LINEAR SOLVER BOOTSTRAP", elapsed1/numRepeats)
    testCases.print("LINEAR SWAP BOOTSTRAP", elapsed2/numRepeats)

###############################################################################


def test_bloombergPricingExample(interp_type):

    """ This is an example of a replication of a BBG example from
    https://github.com/vilen22/curve-building/blob/master/Bloomberg%20Curve%20Building%20Replication.xlsx

    """
    valuation_date = Date(6, 6, 2018)

    # We do the O/N rate which settles on trade date
    spot_days = 0
    settlement_date = valuation_date.addWeekDays(spot_days)
    depoDCCType = DayCountTypes.ACT_360
    depos = []
    deposit_rate = 0.0231381
    maturity_date = settlement_date.addMonths(3)
    depo = IborDeposit(settlement_date, maturity_date, deposit_rate,
                       depoDCCType)
    depos.append(depo)

    futs = []
    fut = IborFuture(valuation_date, 1); futs.append(fut)
    fut = IborFuture(valuation_date, 2); futs.append(fut)
    fut = IborFuture(valuation_date, 3); futs.append(fut)
    fut = IborFuture(valuation_date, 4); futs.append(fut)
    fut = IborFuture(valuation_date, 5); futs.append(fut)
    fut = IborFuture(valuation_date, 6); futs.append(fut)

    fras = [None]*6
    fras[0] = futs[0].toFRA(97.6675, -0.00005)
    fras[1] = futs[1].toFRA(97.5200, -0.00060)
    fras[2] = futs[2].toFRA(97.3550, -0.00146)
    fras[3] = futs[3].toFRA(97.2450, -0.00263)
    fras[4] = futs[4].toFRA(97.1450, -0.00411)
    fras[5] = futs[5].toFRA(97.0750, -0.00589)

    accrual = DayCountTypes.THIRTY_E_360
    freq = FrequencyTypes.SEMI_ANNUAL

    spot_days = 2
    settlement_date = valuation_date.addWeekDays(spot_days)
    notional = ONE_MILLION
    fixed_leg_type = FinSwapTypes.PAY

    swaps = []
    swap = FinIborSwapOLD(settlement_date, "2Y", fixed_leg_type, (2.77417+2.77844)/200, freq, accrual); swaps.append(swap)
    swap = FinIborSwapOLD(settlement_date, "3Y", fixed_leg_type, (2.86098+2.86582)/200, freq, accrual); swaps.append(swap)
    swap = FinIborSwapOLD(settlement_date, "4Y", fixed_leg_type, (2.90240+2.90620)/200, freq, accrual); swaps.append(swap)
    swap = FinIborSwapOLD(settlement_date, "5Y", fixed_leg_type, (2.92944+2.92906)/200, freq, accrual); swaps.append(swap)
    swap = FinIborSwapOLD(settlement_date, "6Y", fixed_leg_type, (2.94001+2.94499)/200, freq, accrual); swaps.append(swap)
    swap = FinIborSwapOLD(settlement_date, "7Y", fixed_leg_type, (2.95352+2.95998)/200, freq, accrual); swaps.append(swap)
    swap = FinIborSwapOLD(settlement_date, "8Y", fixed_leg_type, (2.96830+2.97400)/200, freq, accrual); swaps.append(swap)
    swap = FinIborSwapOLD(settlement_date, "9Y", fixed_leg_type, (2.98403+2.98817)/200, freq, accrual); swaps.append(swap)
    swap = FinIborSwapOLD(settlement_date, "10Y", fixed_leg_type, (2.99716+3.00394)/200, freq, accrual); swaps.append(swap)
    swap = FinIborSwapOLD(settlement_date, "11Y", fixed_leg_type, (3.01344+3.01596)/200, freq, accrual); swaps.append(swap)
    swap = FinIborSwapOLD(settlement_date, "12Y", fixed_leg_type, (3.02276+3.02684)/200, freq, accrual); swaps.append(swap)
    swap = FinIborSwapOLD(settlement_date, "15Y", fixed_leg_type, (3.04092+3.04508)/200, freq, accrual); swaps.append(swap)
    swap = FinIborSwapOLD(settlement_date, "20Y", fixed_leg_type, (3.04417+3.05183)/200, freq, accrual); swaps.append(swap)
    swap = FinIborSwapOLD(settlement_date, "25Y", fixed_leg_type, (3.03219+3.03621)/200, freq, accrual); swaps.append(swap)
    swap = FinIborSwapOLD(settlement_date, "30Y", fixed_leg_type, (3.01030+3.01370)/200, freq, accrual); swaps.append(swap)
    swap = FinIborSwapOLD(settlement_date, "40Y", fixed_leg_type, (2.96946+2.97354)/200, freq, accrual); swaps.append(swap)
    swap = FinIborSwapOLD(settlement_date, "50Y", fixed_leg_type, (2.91552+2.93748)/200, freq, accrual); swaps.append(swap)

    libor_curve = FinIborSingleCurveOLD(valuation_date, depos, fras, swaps, interp_type)

    # The valuation of 53714.55 is very close to the spreadsheet value 53713.96
    principal = 0.0

    testCases.header("VALUATION TO TODAY DATE"," PV")
    testCases.print("VALUE:", swaps[0].value(valuation_date, libor_curve, libor_curve, None))
    testCases.print("FIXED:", swaps[0].fixed_legValue(valuation_date, libor_curve))
    testCases.print("FLOAT:", swaps[0].floatLegValue(valuation_date, libor_curve, libor_curve, None))

    testCases.header("VALUATION TO SWAP SETTLEMENT DATE"," PV")
    testCases.print("VALUE:", swaps[0].value(settlement_date, libor_curve, libor_curve, None))
    testCases.print("FIXED:", swaps[0].fixed_legValue(settlement_date, libor_curve))
    testCases.print("FLOAT:", swaps[0].floatLegValue(settlement_date, libor_curve, libor_curve, None))

    # swaps[0].print_fixed_leg_pv()
    # swaps[0].print_float_leg_pv()

    if 1==0:
        plt.figure()
    
        years = np.linspace(0, 50, 500)    
        dates = settlement_date.addYears(years)
        fwds = libor_curve.fwd(dates)
        plt.plot(years, fwds, label = "Fwd Rate")
        plt.title(interp_type)
        plt.xlabel("Years")
        plt.legend()
    
        years = np.linspace(0, 50, 500)    
        dates = settlement_date.addYears(years)
        fwds = libor_curve.zero_rate(dates)
        plt.plot(years, fwds, label = "Zero Rate")
        plt.title(interp_type)
        plt.xlabel("Years")
        plt.ylabel("Rate")
        plt.legend()
    
###############################################################################

if 1==0:
    for interp_type in FinInterpTypes:
        start = time.time()
        test_bloombergPricingExample(interp_type)
        end = time.time()
        print(interp_type, end - start)

test_bloombergPricingExample(InterpTypes.FLAT_FWD_RATES)
test_derivativePricingExample()
test_FinIborDepositsOnly()
test_FinIborFRAsOnly()
test_FinIborDepositsFRAsSwaps()
test_FinIborDepositsFuturesSwaps()

testCases.compareTestCases()
