###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes

from financepy.products.rates.ibor_fra import IborFRA
from financepy.products.rates.ibor_future import IborFuture
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_swap import IborSwap
from financepy.market.curves.interpolator import InterpTypes

from financepy.utils.math import ONE_MILLION
from financepy.utils.global_types import SwapTypes

from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.dual_curve import IborDualCurve
from financepy.products.rates.ois_curve import OISCurve
from financepy.products.rates.ois import OIS

import numpy as np
import matplotlib.pyplot as plt

test_cases = FinTestCases(__file__, globalTestCaseMode)

PLOT_GRAPHS = False

###############################################################################


def buildOIS(value_dt):
    """ Build the OIS funding curve from futures (FRAs) and OIS """

    spot_days = 0
    spot_days = 0
    settle_dt = value_dt.add_weekdays(spot_days)
    fixed_leg_type = SwapTypes.PAY

    fras = []
    # 1 x 4 FRA

    swaps = []
    fixed_freq_type = FrequencyTypes.SEMI_ANNUAL
    fixed_dcc_type = DayCountTypes.ACT_365F

    swap_rate = 0.000022
    maturity_dt = settle_dt.add_months(24)
    swap = OIS(settle_dt, maturity_dt, fixed_leg_type, swap_rate,
               fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)

    swap_rate += 0.000
    fixed_leg_type = SwapTypes.PAY
    maturity_dt = settle_dt.add_months(36)
    swap = OIS(settle_dt, maturity_dt, fixed_leg_type, swap_rate,
               fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)

    swap_rate += 0.000
    maturity_dt = settle_dt.add_months(48)
    swap = OIS(settle_dt, maturity_dt, fixed_leg_type, swap_rate,
               fixed_freq_type,
               fixed_dcc_type)
    swaps.append(swap)

    swap_rate = 0.02
    maturity_dt = settle_dt.add_months(60)
    swap = OIS(settle_dt, maturity_dt, fixed_leg_type, swap_rate,
               fixed_freq_type,
               fixed_dcc_type)
    swaps.append(swap)

    maturity_dt = settle_dt.add_months(72)
    swap = OIS(settle_dt, maturity_dt, fixed_leg_type, swap_rate,
               fixed_freq_type,
               fixed_dcc_type)
    swaps.append(swap)

    maturity_dt = settle_dt.add_months(84)
    swap = OIS(settle_dt, maturity_dt, fixed_leg_type, swap_rate,
               fixed_freq_type,
               fixed_dcc_type)
    swaps.append(swap)

    maturity_dt = settle_dt.add_months(96)
    swap = OIS(settle_dt, maturity_dt, fixed_leg_type, swap_rate,
               fixed_freq_type,
               fixed_dcc_type)
    swaps.append(swap)

    maturity_dt = settle_dt.add_months(108)
    swap = OIS(settle_dt, maturity_dt, fixed_leg_type, swap_rate,
               fixed_freq_type,
               fixed_dcc_type)
    swaps.append(swap)

    maturity_dt = settle_dt.add_months(120)
    swap = OIS(settle_dt, maturity_dt, fixed_leg_type, swap_rate,
               fixed_freq_type,
               fixed_dcc_type)
    swaps.append(swap)

    maturity_dt = settle_dt.add_months(132)
    swap = OIS(settle_dt, maturity_dt, fixed_leg_type, swap_rate,
               fixed_freq_type,
               fixed_dcc_type)
    swaps.append(swap)

    maturity_dt = settle_dt.add_months(144)
    swap = OIS(settle_dt, maturity_dt, fixed_leg_type, swap_rate,
               fixed_freq_type,
               fixed_dcc_type)
    swaps.append(swap)

    maturity_dt = settle_dt.add_months(180)
    swap = OIS(settle_dt, maturity_dt, fixed_leg_type, swap_rate,
               fixed_freq_type,
               fixed_dcc_type)
    swaps.append(swap)

    maturity_dt = settle_dt.add_months(240)
    swap = OIS(settle_dt, maturity_dt, fixed_leg_type, swap_rate,
               fixed_freq_type,
               fixed_dcc_type)
    swaps.append(swap)

    maturity_dt = settle_dt.add_months(300)
    swap = OIS(settle_dt, maturity_dt, fixed_leg_type, swap_rate,
               fixed_freq_type,
               fixed_dcc_type)
    swaps.append(swap)

    maturity_dt = settle_dt.add_months(360)
    swap = OIS(settle_dt, maturity_dt, fixed_leg_type, swap_rate,
               fixed_freq_type,
               fixed_dcc_type)
    swaps.append(swap)

    oisCurve = OISCurve(value_dt,
                        [],
                        fras,
                        swaps)

    return oisCurve

###############################################################################


def test_bloombergPricingExample():
    """ This is an example of a replication of a BBG example from
    https://github.com/vilen22/curve-building/blob/master/Bloomberg%20Curve%20Building%20Replication.xlsx

    """
    value_dt = Date(6, 6, 2018)

    # We do the O/N rate which settles on trade date
    spot_days = 0
    settle_dt = value_dt.add_weekdays(spot_days)
    depoDCCType = DayCountTypes.ACT_360
    depos = []
    deposit_rate = 0.0231381
    maturity_dt = settle_dt.add_months(3)
    depo = IborDeposit(settle_dt, maturity_dt, deposit_rate,
                       depoDCCType)
    depos.append(depo)

    futs = []
    fut = IborFuture(value_dt, 1)
    futs.append(fut)
    fut = IborFuture(value_dt, 2)
    futs.append(fut)
    fut = IborFuture(value_dt, 3)
    futs.append(fut)
    fut = IborFuture(value_dt, 4)
    futs.append(fut)
    fut = IborFuture(value_dt, 5)
    futs.append(fut)
    fut = IborFuture(value_dt, 6)
    futs.append(fut)

    fras = [None]*6
    fras[0] = futs[0].to_fra(97.6675, -0.00005)
    fras[1] = futs[1].to_fra(97.5200, -0.00060)
    fras[2] = futs[2].to_fra(97.3550, -0.00146)
    fras[3] = futs[3].to_fra(97.2450, -0.00263)
    fras[4] = futs[4].to_fra(97.1450, -0.00411)
    fras[5] = futs[5].to_fra(97.0750, -0.00589)

    accrual = DayCountTypes.THIRTY_E_360
    freq = FrequencyTypes.SEMI_ANNUAL

    spot_days = 2
    settle_dt = value_dt.add_weekdays(spot_days)
    fixed_leg_type = SwapTypes.PAY
    interp_type = InterpTypes.FLAT_FWD_RATES

    swaps = []
    swap = IborSwap(settle_dt, "2Y", fixed_leg_type,
                    (2.77417 + 2.77844) / 200, freq, accrual)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "3Y", fixed_leg_type,
                    (2.86098 + 2.86582) / 200, freq, accrual)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "4Y", fixed_leg_type,
                    (2.90240 + 2.90620) / 200, freq, accrual)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "5Y", fixed_leg_type,
                    (2.92944 + 2.92906) / 200, freq, accrual)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "6Y", fixed_leg_type,
                    (2.94001 + 2.94499) / 200, freq, accrual)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "7Y", fixed_leg_type,
                    (2.95352 + 2.95998) / 200, freq, accrual)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "8Y", fixed_leg_type,
                    (2.96830 + 2.97400) / 200, freq, accrual)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "9Y", fixed_leg_type,
                    (2.98403 + 2.98817) / 200, freq, accrual)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "10Y", fixed_leg_type,
                    (2.99716 + 3.00394) / 200, freq, accrual)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "11Y", fixed_leg_type,
                    (3.01344 + 3.01596) / 200, freq, accrual)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "12Y", fixed_leg_type,
                    (3.02276 + 3.02684) / 200, freq, accrual)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "15Y", fixed_leg_type,
                    (3.04092 + 3.04508) / 200, freq, accrual)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "20Y", fixed_leg_type,
                    (3.04417 + 3.05183) / 200, freq, accrual)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "25Y", fixed_leg_type,
                    (3.03219 + 3.03621) / 200, freq, accrual)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "30Y", fixed_leg_type,
                    (3.01030 + 3.01370) / 200, freq, accrual)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "40Y", fixed_leg_type,
                    (2.96946 + 2.97354) / 200, freq, accrual)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "50Y", fixed_leg_type,
                    (2.91552 + 2.93748) / 200, freq, accrual)
    swaps.append(swap)

    libor_curve = IborSingleCurve(
        value_dt, depos, fras, swaps, interp_type, True)

    test_cases.banner("======================================================")
    test_cases.banner("SINGLE CURVE VALUATION")
    test_cases.header("LABEL", "VALUE")
    test_cases.print("VALUE:", swaps[0].value(
        value_dt, libor_curve, libor_curve, None))
    test_cases.print("FIXED:", swaps[0].fixed_leg.value(
        value_dt, libor_curve))
    test_cases.print("FLOAT:", swaps[0].float_leg.value(
        value_dt, libor_curve, libor_curve, None))

    test_cases.banner("======================================================")
    test_cases.banner("SINGLE CURVE VALUATION TO SWAP SETTLEMENT DATE")
    test_cases.header("LABEL", "VALUE")
    test_cases.print("VALUE:", swaps[0].value(
        settle_dt, libor_curve, libor_curve, None))
    test_cases.print("FIXED:", swaps[0].fixed_leg.value(
        settle_dt, libor_curve))
    test_cases.print("FLOAT:", swaps[0].float_leg.value(
        settle_dt, libor_curve, libor_curve, None))
    test_cases.banner("======================================================")

#    swaps[0].print_fixed_leg_pv()
#    swaps[0].print_float_leg_pv()

    oisCurve = buildOIS(value_dt)
#    print(oisCurve)

    liborDualCurve = IborDualCurve(value_dt, oisCurve, depos, fras, swaps,
                                   InterpTypes.FLAT_FWD_RATES, True)
#    print(liborDualCurve)

    # The valuation of 53714.55 is very close to the spreadsheet value 53713.96

    test_cases.header("VALUATION TO TODAY DATE", " PV")
    test_cases.print("VALUE:", swaps[0].value(
        value_dt, oisCurve, liborDualCurve, None))
    test_cases.print("FIXED:", swaps[0].fixed_leg.value(
        value_dt, oisCurve))
    test_cases.print("FLOAT:", swaps[0].float_leg.value(
        value_dt, oisCurve, libor_curve, None))

    test_cases.header("VALUATION TO SWAP SETTLEMENT DATE", " PV")
    test_cases.print("VALUE:", swaps[0].value(
        settle_dt, oisCurve, liborDualCurve, None))
    test_cases.print("FIXED:", swaps[0].fixed_leg.value(
        settle_dt, oisCurve))
    test_cases.print("FLOAT:", swaps[0].float_leg.value(
        settle_dt, oisCurve, liborDualCurve, None, ))

#    swaps[0].print_fixed_leg_pv()
#    swaps[0].print_float_leg_pv()

    PLOT = False
    if PLOT is True:

        years = np.linspace(0, 5, 21)
        dates = settle_dt.add_years(years)

        singleCurveFwds = libor_curve.fwd(dates)
        plt.plot(years, singleCurveFwds, label="Single Libor Curve")

        oisCurveFwds = oisCurve.fwd(dates)
        plt.plot(years, oisCurveFwds, label="OIS Curve")

        index_curveFwds = liborDualCurve.fwd(dates)
        plt.plot(years, index_curveFwds, label="Libor Index Curve")

        plt.legend()

    # swaps[0].print_fixed_leg_pv()
    # swaps[0].print_float_leg_pv()

###############################################################################


def test_swapValuationExample():

    # Example from
    # https://blog.deriscope.com/index.php/en/excel-interest-rate-swap-price-dual-bootstrapping-curve

    vBloomberg = 388147

    value_dt = Date(30, 11, 2018)

    start_dt = Date(27, 12, 2017)
    maturity_dt = Date(27, 12, 2067)
    notional = 10 * ONE_MILLION
    fixed_leg_type = SwapTypes.RECEIVE

    fixedRate = 0.0150
    fixed_dcc_type = DayCountTypes.THIRTY_360_BOND
    fixed_freq_type = FrequencyTypes.ANNUAL

    float_spread = 0.0
    float_dcc_type = DayCountTypes.ACT_360
    float_freq_type = FrequencyTypes.SEMI_ANNUAL

    offMarketSwap = IborSwap(start_dt, maturity_dt, fixed_leg_type,
                             fixedRate, fixed_freq_type, fixed_dcc_type,
                             notional,
                             float_spread, float_freq_type, float_dcc_type)

    interp_type = InterpTypes.LINEAR_ZERO_RATES

    depoDCCType = DayCountTypes.ACT_360
    depos = []

    ###########################################################################
    # MARKET
    ###########################################################################

    spot_days = 0
    settle_dt = value_dt.add_weekdays(spot_days)
    depo = IborDeposit(settle_dt, "6M", -0.2510 / 100.0, depoDCCType)
    depos.append(depo)

    fras = []
    fraDCCType = DayCountTypes.ACT_360

    fra = IborFRA(settle_dt.add_tenor("1M"),
                  "6M", -0.2450 / 100.0, fraDCCType)
    fras.append(fra)
    fra = IborFRA(settle_dt.add_tenor("2M"),
                  "6M", -0.2435 / 100.0, fraDCCType)
    fras.append(fra)
    fra = IborFRA(settle_dt.add_tenor("3M"),
                  "6M", -0.2400 / 100.0, fraDCCType)
    fras.append(fra)
    fra = IborFRA(settle_dt.add_tenor("4M"),
                  "6M", -0.2360 / 100.0, fraDCCType)
    fras.append(fra)
    fra = IborFRA(settle_dt.add_tenor("5M"),
                  "6M", -0.2285 / 100.0, fraDCCType)
    fras.append(fra)
    fra = IborFRA(settle_dt.add_tenor("6M"),
                  "6M", -0.2230 / 100.0, fraDCCType)
    fras.append(fra)
    fra = IborFRA(settle_dt.add_tenor("7M"),
                  "6M", -0.2110 / 100.0, fraDCCType)
    fras.append(fra)
    fra = IborFRA(settle_dt.add_tenor("8M"),
                  "6M", -0.1990 / 100.0, fraDCCType)
    fras.append(fra)
    fra = IborFRA(settle_dt.add_tenor("9M"),
                  "6M", -0.1850 / 100.0, fraDCCType)
    fras.append(fra)
    fra = IborFRA(settle_dt.add_tenor("10M"),
                  "6M", -0.1680 / 100.0, fraDCCType)
    fras.append(fra)
    fra = IborFRA(settle_dt.add_tenor("11M"),
                  "6M", -0.1510 / 100.0, fraDCCType)
    fras.append(fra)
    fra = IborFRA(settle_dt.add_tenor("12M"),
                  "6M", -0.1360 / 100.0, fraDCCType)
    fras.append(fra)

    swaps = []
    fixed_leg_type = SwapTypes.PAY
    fixed_dcc_type = DayCountTypes.THIRTY_360_BOND
    fixed_freq_type = FrequencyTypes.ANNUAL

    swap = IborSwap(settle_dt, "2Y", fixed_leg_type, -
                    0.1525 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "3Y", fixed_leg_type, -
                    0.0185 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "4Y", fixed_leg_type,
                    0.1315 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "5Y", fixed_leg_type,
                    0.2745 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "6Y", fixed_leg_type,
                    0.4135 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "7Y", fixed_leg_type,
                    0.5439 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "8Y", fixed_leg_type,
                    0.6652 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "9Y", fixed_leg_type,
                    0.7784 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "10Y", fixed_leg_type,
                    0.8799 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "11Y", fixed_leg_type,
                    0.9715 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "12Y", fixed_leg_type,
                    1.0517 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "15Y", fixed_leg_type,
                    1.2369 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "20Y", fixed_leg_type,
                    1.3965 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "25Y", fixed_leg_type,
                    1.4472 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "30Y", fixed_leg_type,
                    1.4585 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "35Y", fixed_leg_type,
                    1.4595 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "40Y", fixed_leg_type,
                    1.4535 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "45Y", fixed_leg_type,
                    1.4410 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = IborSwap(settle_dt, "50Y", fixed_leg_type,
                    1.4335 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)

    iborDepos = depos.copy()
    iborFras = fras.copy()
    ibor_swaps = swaps.copy()

    iborCurve = IborSingleCurve(
        value_dt, iborDepos, iborFras, ibor_swaps, interp_type)
    v1 = offMarketSwap.value(value_dt, iborCurve,
                             iborCurve, -0.268/100.0)

    test_cases.banner("DERISCOPE EXAMPLE REPLICATION")
    test_cases.header("LABEL", "VALUE")
    test_cases.print("BBG VALUE", vBloomberg)
    test_cases.print("FP ONE CURVE VALUE", v1)

    ###############################################################################

    depoDCCType = DayCountTypes.ACT_360
    depos = []

    spot_days = 0
    settle_dt = value_dt.add_weekdays(spot_days)
    depo = IborDeposit(settle_dt, "1D", -0.3490 / 100.0, depoDCCType)
    depos.append(depo)

    fras = []

    swaps = []
    fixed_leg_type = SwapTypes.PAY
    fixed_dcc_type = DayCountTypes.ACT_365F
    fixed_freq_type = FrequencyTypes.ANNUAL

    # Standard OIS with standard annual terms
    swap = OIS(settle_dt, "2W", fixed_leg_type, -
               0.3600 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "1M", fixed_leg_type, -
               0.3560 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "2M", fixed_leg_type, -
               0.3570 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "3M", fixed_leg_type, -
               0.3580 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "4M", fixed_leg_type, -
               0.3575 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "5M", fixed_leg_type, -
               0.3578 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "6M", fixed_leg_type, -
               0.3580 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "7M", fixed_leg_type, -
               0.3600 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "8M", fixed_leg_type, -
               0.3575 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "9M", fixed_leg_type, -
               0.3569 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "10M", fixed_leg_type, -
               0.3553 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "11M", fixed_leg_type, -
               0.3534 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "12M", fixed_leg_type, -
               0.3496 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "18M", fixed_leg_type, -
               0.3173 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)

    swap = OIS(settle_dt, "2Y", fixed_leg_type, -
               0.2671 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "30M", fixed_leg_type, -
               0.2070 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "3Y", fixed_leg_type, -
               0.1410 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "4Y", fixed_leg_type, -
               0.0060 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "5Y", fixed_leg_type,
               0.1285 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "6Y", fixed_leg_type,
               0.2590 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "7Y", fixed_leg_type,
               0.3830 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "8Y", fixed_leg_type,
               0.5020 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "9Y", fixed_leg_type,
               0.6140 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "10Y", fixed_leg_type,
               0.7160 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "11Y", fixed_leg_type,
               0.8070 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "12Y", fixed_leg_type,
               0.8890 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "15Y", fixed_leg_type,
               1.0790 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "20Y", fixed_leg_type,
               1.2460 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "25Y", fixed_leg_type,
               1.3055 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "30Y", fixed_leg_type,
               1.3270 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "35Y", fixed_leg_type,
               1.3315 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "40Y", fixed_leg_type,
               1.3300 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)
    swap = OIS(settle_dt, "50Y", fixed_leg_type,
               1.3270 / 100.0, fixed_freq_type, fixed_dcc_type)
    swaps.append(swap)

    oisDepos = depos.copy()
    oisFras = fras.copy()
    oisSwaps = swaps.copy()

#    oisCurveFF = OISCurve(value_dt, oisDepos, oisFras, oisSwaps, interp_type)

    iborDualCurve = IborDualCurve(
        value_dt, oisCurveFF, iborDepos, iborFras, ibor_swaps, interp_type)

#    v2 = offMarketSwap.value(value_dt, oisCurveFF, iborDualCurve, -0.268/100.0)

#    test_cases.print("FP DUAL CURVE VALUE", v2)

#    swap_rate = offMarketSwap.swap_rate(value_dt, oisCurveFF, iborCurve, -0.268/100.0)

#    test_cases.print("FP DUAL CURVE SWAP RATE", swap_rate)

#    offMarketSwap.printFixedLegFlows()
#    offMarketSwap.printFloatLegFlows()
#    offMarketSwap.print_fixed_leg_pv()
#    offMarketSwap.print_float_leg_pv()


###############################################################################

# test_swapValuationExample()

test_bloombergPricingExample()

test_cases.compareTestCases()
