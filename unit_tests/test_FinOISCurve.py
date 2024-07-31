###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.global_types import SwapTypes
from financepy.market.curves.interpolator import InterpTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ois_curve import OISCurve
from financepy.products.rates.ois import OIS
from financepy.products.rates.ibor_future import IborFuture
from financepy.products.rates.ibor_fra import IborFRA
from financepy.utils.calendar import CalendarTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
import matplotlib.pyplot as plt
import numpy as np


def test_bloombergPricingExample():
    """ This is an example of a replication of a BBG example from
    https://github.com/vilen22/curve-building/blob/master/Bloomberg%20Curve%20Building%20Replication.xlsx
    """

    value_dt = Date(6, 6, 2018)

    # We do the O/N rate which settles on trade date
    spot_days = 0
    settle_dt = value_dt.add_weekdays(spot_days)
    accrual = DayCountTypes.THIRTY_E_360

    depo = IborDeposit(settle_dt, "1D", 1.712 / 100.0, accrual)
    depos = [depo]

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
    payRec = SwapTypes.PAY
    lag = 1  # Not used

    swaps = []
    swap = OIS(settle_dt, "2Y", payRec,
               (2.77417 + 2.77844) / 200, freq, accrual)
    swaps.append(swap)
    swap = OIS(settle_dt, "3Y", payRec,
               (2.86098 + 2.86582) / 200, freq, accrual)
    swaps.append(swap)
    swap = OIS(settle_dt, "4Y", payRec,
               (2.90240 + 2.90620) / 200, freq, accrual)
    swaps.append(swap)
    swap = OIS(settle_dt, "5Y", payRec,
               (2.92944 + 2.92906) / 200, freq, accrual)
    swaps.append(swap)
    swap = OIS(settle_dt, "6Y", payRec,
               (2.94001 + 2.94499) / 200, freq, accrual)
    swaps.append(swap)
    swap = OIS(settle_dt, "7Y", payRec,
               (2.95352 + 2.95998) / 200, freq, accrual)
    swaps.append(swap)
    swap = OIS(settle_dt, "8Y", payRec,
               (2.96830 + 2.97400) / 200, freq, accrual)
    swaps.append(swap)
    swap = OIS(settle_dt, "9Y", payRec,
               (2.98403 + 2.98817) / 200, freq, accrual)
    swaps.append(swap)
    swap = OIS(settle_dt, "10Y", payRec,
               (2.99716 + 3.00394) / 200, freq, accrual)
    swaps.append(swap)
    swap = OIS(settle_dt, "11Y", payRec,
               (3.01344 + 3.01596) / 200, freq, accrual)
    swaps.append(swap)
    swap = OIS(settle_dt, "12Y", payRec,
               (3.02276 + 3.02684) / 200, freq, accrual)
    swaps.append(swap)
    swap = OIS(settle_dt, "15Y", payRec,
               (3.04092 + 3.04508) / 200, freq, accrual)
    swaps.append(swap)
    swap = OIS(settle_dt, "20Y", payRec,
               (3.04417 + 3.05183) / 200, freq, accrual)
    swaps.append(swap)
    swap = OIS(settle_dt, "25Y", payRec,
               (3.03219 + 3.03621) / 200, freq, accrual)
    swaps.append(swap)
    swap = OIS(settle_dt, "30Y", payRec,
               (3.01030 + 3.01370) / 200, freq, accrual)
    swaps.append(swap)
    swap = OIS(settle_dt, "40Y", payRec,
               (2.96946 + 2.97354) / 200, freq, accrual)
    swaps.append(swap)
    swap = OIS(settle_dt, "50Y", payRec,
               (2.91552 + 2.93748) / 200, freq, accrual)
    swaps.append(swap)

    oisCurve = OISCurve(value_dt, depos, fras, swaps)

    # The valuation of 53714.55 is very close to the spreadsheet value 53713.96
    assert round(swaps[0].value(value_dt, oisCurve, None), 4) == 0.0
    assert round(-swaps[0].fixed_leg.value(value_dt,
                                            oisCurve), 4) == 53708.2780
    assert round(swaps[0].float_leg.value(
        value_dt, oisCurve, None), 4) == 53708.2780

    assert round(swaps[0].value(settle_dt, oisCurve, None), 4) == 0.0
    assert round(-swaps[0].fixed_leg.value(settle_dt,
                                            oisCurve), 4) == 53714.3020
    assert round(swaps[0].float_leg.value(
        settle_dt, oisCurve, None), 4) == 53714.3020
