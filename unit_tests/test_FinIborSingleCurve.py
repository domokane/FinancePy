###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################
import pytest

from financepy.utils.global_types import SwapTypes
from financepy.utils.math import ONE_MILLION
from financepy.market.curves.interpolator import InterpTypes
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.rates.ibor_fra import IborFRA
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_future import IborFuture
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
from financepy.utils.calendar import Calendar, CalendarTypes


def test_bloombergPricingExample():
    """This is an example of a replication of a BBG example from
    https://github.com/vilen22/curve-building/blob/master/Bloomberg%20Curve%20Building%20Replication.xlsx

    """
    value_dt = Date(6, 6, 2018)
    interp_type = InterpTypes.FLAT_FWD_RATES

    # We do the O/N rate which settles on trade date
    spot_days = 0
    settle_dt = value_dt.add_weekdays(spot_days)
    depoDCCType = DayCountTypes.ACT_360
    depos = []
    deposit_rate = 0.0231381
    maturity_dt = settle_dt.add_months(3)
    depo = IborDeposit(settle_dt, maturity_dt, deposit_rate, depoDCCType)
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

    fras = [None] * 6
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
    notional = ONE_MILLION
    fixed_leg_type = SwapTypes.PAY

    swaps = []
    swap = IborSwap(
        settle_dt,
        "2Y",
        fixed_leg_type,
        (2.77417 + 2.77844) / 200,
        freq,
        accrual,
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt,
        "3Y",
        fixed_leg_type,
        (2.86098 + 2.86582) / 200,
        freq,
        accrual,
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt,
        "4Y",
        fixed_leg_type,
        (2.90240 + 2.90620) / 200,
        freq,
        accrual,
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt,
        "5Y",
        fixed_leg_type,
        (2.92944 + 2.92906) / 200,
        freq,
        accrual,
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt,
        "6Y",
        fixed_leg_type,
        (2.94001 + 2.94499) / 200,
        freq,
        accrual,
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt,
        "7Y",
        fixed_leg_type,
        (2.95352 + 2.95998) / 200,
        freq,
        accrual,
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt,
        "8Y",
        fixed_leg_type,
        (2.96830 + 2.97400) / 200,
        freq,
        accrual,
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt,
        "9Y",
        fixed_leg_type,
        (2.98403 + 2.98817) / 200,
        freq,
        accrual,
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt,
        "10Y",
        fixed_leg_type,
        (2.99716 + 3.00394) / 200,
        freq,
        accrual,
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt,
        "11Y",
        fixed_leg_type,
        (3.01344 + 3.01596) / 200,
        freq,
        accrual,
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt,
        "12Y",
        fixed_leg_type,
        (3.02276 + 3.02684) / 200,
        freq,
        accrual,
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt,
        "15Y",
        fixed_leg_type,
        (3.04092 + 3.04508) / 200,
        freq,
        accrual,
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt,
        "20Y",
        fixed_leg_type,
        (3.04417 + 3.05183) / 200,
        freq,
        accrual,
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt,
        "25Y",
        fixed_leg_type,
        (3.03219 + 3.03621) / 200,
        freq,
        accrual,
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt,
        "30Y",
        fixed_leg_type,
        (3.01030 + 3.01370) / 200,
        freq,
        accrual,
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt,
        "40Y",
        fixed_leg_type,
        (2.96946 + 2.97354) / 200,
        freq,
        accrual,
    )
    swaps.append(swap)
    swap = IborSwap(
        settle_dt,
        "50Y",
        fixed_leg_type,
        (2.91552 + 2.93748) / 200,
        freq,
        accrual,
    )
    swaps.append(swap)

    libor_curve = IborSingleCurve(value_dt, depos, fras, swaps, interp_type)

    # The valuation of 53714.55 is very close to the spreadsheet value 53713.96
    principal = 0.0

    # Pay fixed so make fixed leg value negative
    assert (
        round(swaps[0].value(value_dt, libor_curve, libor_curve, None), 4)
        == 0.0
    )
    assert (
        round(-swaps[0].fixed_leg.value(value_dt, libor_curve), 4)
        == 53707.6667
    )
    assert (
        round(
            swaps[0].float_leg.value(value_dt, libor_curve, libor_curve, None),
            4,
        )
        == 53707.6667
    )

    # Pay fixed so make fixed leg value negative
    assert (
        round(swaps[0].value(settle_dt, libor_curve, libor_curve, None), 4)
        == 0.0
    )
    assert (
        round(-swaps[0].fixed_leg.value(settle_dt, libor_curve), 4)
        == 53714.5507
    )
    assert (
        round(
            swaps[0].float_leg.value(
                settle_dt, libor_curve, libor_curve, None
            ),
            4,
        )
        == 53714.5507
    )


@pytest.mark.parametrize("interp_type", InterpTypes)
def test_RepriceInputsForAllInterpChoices(interp_type):
    valuation_date = Date(6, 10, 2001)

    cal = CalendarTypes.UNITED_KINGDOM

    depoDCCType = DayCountTypes.ACT_360
    depos = []
    spot_days = 2
    settlement_date = valuation_date.add_weekdays(spot_days)
    depo = IborDeposit(
        settlement_date, "3M", 4.2 / 100.0, depoDCCType, cal_type=cal
    )
    depos.append(depo)

    fraDCCType = DayCountTypes.ACT_360
    fras = []
    fra = IborFRA(
        settlement_date.add_tenor("3M"),
        "3M",
        4.20 / 100.0,
        fraDCCType,
        cal_type=cal,
    )
    fras.append(fra)

    swaps = []
    swapType = SwapTypes.PAY
    fixedDCCType = DayCountTypes.THIRTY_E_360_ISDA
    fixed_freqType = FrequencyTypes.SEMI_ANNUAL

    swap = IborSwap(
        settlement_date,
        "1Y",
        swapType,
        4.20 / 100.0,
        fixed_freqType,
        fixedDCCType,
        cal_type=cal,
    )
    swaps.append(swap)
    swap = IborSwap(
        settlement_date,
        "2Y",
        swapType,
        4.30 / 100.0,
        fixed_freqType,
        fixedDCCType,
        cal_type=cal,
    )
    swaps.append(swap)
    swap = IborSwap(
        settlement_date,
        "3Y",
        swapType,
        4.70 / 100.0,
        fixed_freqType,
        fixedDCCType,
        cal_type=cal,
    )
    swaps.append(swap)
    swap = IborSwap(
        settlement_date,
        "5Y",
        swapType,
        5.40 / 100.0,
        fixed_freqType,
        fixedDCCType,
        cal_type=cal,
    )
    swaps.append(swap)
    swap = IborSwap(
        settlement_date,
        "7Y",
        swapType,
        5.70 / 100.0,
        fixed_freqType,
        fixedDCCType,
        cal_type=cal,
    )
    swaps.append(swap)
    swap = IborSwap(
        settlement_date,
        "10Y",
        swapType,
        6.00 / 100.0,
        fixed_freqType,
        fixedDCCType,
        cal_type=cal,
    )
    swaps.append(swap)
    swap = IborSwap(
        settlement_date,
        "12Y",
        swapType,
        6.10 / 100.0,
        fixed_freqType,
        fixedDCCType,
        cal_type=cal,
    )
    swaps.append(swap)
    swap = IborSwap(
        settlement_date,
        "15Y",
        swapType,
        5.90 / 100.0,
        fixed_freqType,
        fixedDCCType,
        cal_type=cal,
    )
    swaps.append(swap)
    swap = IborSwap(
        settlement_date,
        "20Y",
        swapType,
        5.60 / 100.0,
        fixed_freqType,
        fixedDCCType,
        cal_type=cal,
    )
    swaps.append(swap)
    swap = IborSwap(
        settlement_date,
        "25Y",
        swapType,
        5.55 / 100.0,
        fixed_freqType,
        fixedDCCType,
        cal_type=cal,
    )
    swaps.append(swap)

    optional_interp_params = {
        "sigma": 5.0
    }  # only relevant for interp_type == InterpTypes.TENSION_ZERO_RATES

    # this will throw if inputs are not repriced
    IborSingleCurve(
        valuation_date,
        depos,
        fras,
        swaps,
        interp_type,
        check_refit=True,
        **optional_interp_params
    )

    # If no exception, we are good
    assert True
