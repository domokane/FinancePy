import pytest
import pandas as pd

from financepy.utils.global_types import SwapTypes
from financepy.utils.math import ONE_MILLION
from financepy.market.curves.interpolator import InterpTypes
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.rates.ibor_fra import IborFRA
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_future import IborFuture
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_single_curve_smoothing_calibrator import (
    IborSingleCurveSmoothingCalibrator,
)
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
from financepy.utils.calendar import Calendar, CalendarTypes


@pytest.mark.parametrize("interp_type", InterpTypes)
def test_SmoothFitSimple(interp_type, report_progress=False):
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

    # Create but do not build the initial curve
    do_build = False
    init_curve = IborSingleCurve(
        valuation_date,
        depos,
        fras,
        swaps,
        interp_type,
        check_refit=False,
        do_build=do_build,
    )

    calibrator = IborSingleCurveSmoothingCalibrator(init_curve)

    smooth_param = 1.0
    curve, report = calibrator.fit(
        smoothness=smooth_param, report_progress=report_progress
    )

    if report_progress:
        with pd.option_context(
            "display.max_rows",
            None,
            "display.max_columns",
            None,
            "display.float_format",
            lambda x: f"{x:.4f}",
        ):

            print(report)


@pytest.mark.parametrize("interp_type", [InterpTypes.FLAT_FWD_RATES])
def test_SmoothFit(interp_type, report_progress=False):
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

    # Create but do not build the initial curve
    do_build = False
    init_curve = IborSingleCurve(
        valuation_date,
        depos,
        fras,
        swaps,
        interp_type,
        check_refit=False,
        do_build=do_build,
        **optional_interp_params,
    )

    calibrator = IborSingleCurveSmoothingCalibrator(init_curve)

    # here we go
    smooth_param = 1e-2
    curve, report = calibrator.fit(
        smoothness=smooth_param, report_progress=report_progress
    )
    if report_progress:
        with pd.option_context(
            "display.max_rows",
            None,
            "display.max_columns",
            None,
            "display.float_format",
            lambda x: f"{x:.4f}",
        ):

            print(report)

    # If no exception, we are good


if __name__ == "__main__":
    report_progress = True
    test_SmoothFitSimple(InterpTypes.LINEAR_ZERO_RATES, report_progress)
    test_SmoothFit(InterpTypes.FLAT_FWD_RATES, report_progress)
    # test_SmoothFit(InterpTypes.FINCUBIC_ZERO_RATES, report_progress)
