###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.products.rates.ois import OIS
from financepy.products.rates.ois_curve import OISCurve
from financepy.products.rates.dual_curve import IborDualCurve
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.utils.global_types import SwapTypes
from financepy.utils.math import ONE_MILLION
from financepy.market.curves.interpolator import InterpTypes
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_future import IborFuture
from financepy.products.rates.ibor_fra import IborFRA
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
import matplotlib.pyplot as plt
import numpy as np


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

    assert round(swaps[0].value(
        value_dt, libor_curve, libor_curve, None), 4) == 0.0
    assert round(swaps[0].fixed_leg.value(
        value_dt, libor_curve), 4) == -53707.6667
    assert round(swaps[0].float_leg.value(
        value_dt, libor_curve, libor_curve, None), 4) == 53707.6667

    assert round(swaps[0].value(
        settle_dt, libor_curve, libor_curve, None), 4) == 0.0
    assert round(swaps[0].fixed_leg.value(
        settle_dt, libor_curve), 4) == -53714.5507
    assert round(swaps[0].float_leg.value(
        settle_dt, libor_curve, libor_curve, None), 4) == 53714.5507

    oisCurve = buildOIS(value_dt)

    liborDualCurve = IborDualCurve(value_dt, oisCurve, depos, fras, swaps,
                                   InterpTypes.FLAT_FWD_RATES, True)

    assert round(swaps[0].value(
        value_dt, oisCurve, liborDualCurve, None), 4) == 0.0
    assert round(swaps[0].fixed_leg.value(
        value_dt, oisCurve), 4) == -55524.5642
    assert round(swaps[0].float_leg.value(
        value_dt, oisCurve, liborDualCurve, None), 4) == 55524.5642

    assert round(swaps[0].value(
        settle_dt, oisCurve, liborDualCurve, None), 4) == 0.0
    assert round(swaps[0].fixed_leg.value(
        settle_dt, oisCurve), 4) == -55524.5709
    assert round(swaps[0].float_leg.value(
        settle_dt, oisCurve, liborDualCurve, None), 4) == 55524.5709
