import pandas as pd
import pytest

from financepy.utils.global_types import SwapTypes
from financepy.market.curves.interpolator import InterpTypes
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.rates.ibor_fra import IborFRA
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_single_curve_smoothing_calibrator import (
    IborSingleCurveSmoothingCalibrator,
)
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date

REPORT_PROGRESS = True


@pytest.mark.parametrize("interp_type", InterpTypes)

########################################################################################


def test_smooth_fit_simple(interp_type):

    valuation_date = Date(6, 10, 2001)
    cal = CalendarTypes.UNITED_KINGDOM

    depo_dcc_type = DayCountTypes.ACT_360
    depos = []
    spot_days = 2
    settlement_date = valuation_date.add_weekdays(spot_days)
    depo = IborDeposit(settlement_date, "3M", 4.2 / 100.0, depo_dcc_type, cal_type=cal)
    depos.append(depo)

    fra_dcc_type = DayCountTypes.ACT_360
    fras = []
    fra = IborFRA(
        settlement_date.add_tenor("3M"),
        "3M",
        4.20 / 100.0,
        fra_dcc_type,
        cal_type=cal,
    )
    fras.append(fra)

    swaps = []
    swap_type = SwapTypes.PAY
    fixed_dcc_type = DayCountTypes.THIRTY_E_360_ISDA
    fixed_freq_type = FrequencyTypes.SEMI_ANNUAL

    swap = IborSwap(
        settlement_date,
        "1Y",
        swap_type,
        4.20 / 100.0,
        fixed_freq_type,
        fixed_dcc_type,
        cal_type=cal,
    )
    swaps.append(swap)
    swap = IborSwap(
        settlement_date,
        "2Y",
        swap_type,
        4.30 / 100.0,
        fixed_freq_type,
        fixed_dcc_type,
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
        check_refit_flag=False,
        do_build=do_build,
    )

    calibrator = IborSingleCurveSmoothingCalibrator(init_curve)

    smooth_param = 1.0
    curve, report = calibrator.fit(
        smoothness=smooth_param, report_progress=REPORT_PROGRESS
    )

    if REPORT_PROGRESS:
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

########################################################################################


def test_smooth_fit(interp_type):

    valuation_date = Date(6, 10, 2001)
    cal = CalendarTypes.UNITED_KINGDOM

    depo_dcc_type = DayCountTypes.ACT_360
    depos = []
    spot_days = 2
    settlement_date = valuation_date.add_weekdays(spot_days)
    depo = IborDeposit(settlement_date, "3M", 4.2 / 100.0, depo_dcc_type, cal_type=cal)
    depos.append(depo)

    fra_dcc_type = DayCountTypes.ACT_360
    fras = []
    fra = IborFRA(
        settlement_date.add_tenor("3M"),
        "3M",
        4.20 / 100.0,
        fra_dcc_type,
        cal_type=cal,
    )
    fras.append(fra)

    swaps = []
    swap_type = SwapTypes.PAY
    fixed_dcc_type = DayCountTypes.THIRTY_E_360_ISDA
    fixed_freq_type = FrequencyTypes.SEMI_ANNUAL

    swap = IborSwap(
        settlement_date,
        "1Y",
        swap_type,
        4.20 / 100.0,
        fixed_freq_type,
        fixed_dcc_type,
        cal_type=cal,
    )
    swaps.append(swap)
    swap = IborSwap(
        settlement_date,
        "2Y",
        swap_type,
        4.30 / 100.0,
        fixed_freq_type,
        fixed_dcc_type,
        cal_type=cal,
    )
    swaps.append(swap)
    swap = IborSwap(
        settlement_date,
        "3Y",
        swap_type,
        4.70 / 100.0,
        fixed_freq_type,
        fixed_dcc_type,
        cal_type=cal,
    )
    swaps.append(swap)
    swap = IborSwap(
        settlement_date,
        "5Y",
        swap_type,
        5.40 / 100.0,
        fixed_freq_type,
        fixed_dcc_type,
        cal_type=cal,
    )
    swaps.append(swap)
    swap = IborSwap(
        settlement_date,
        "7Y",
        swap_type,
        5.70 / 100.0,
        fixed_freq_type,
        fixed_dcc_type,
        cal_type=cal,
    )
    swaps.append(swap)
    swap = IborSwap(
        settlement_date,
        "10Y",
        swap_type,
        6.00 / 100.0,
        fixed_freq_type,
        fixed_dcc_type,
        cal_type=cal,
    )
    swaps.append(swap)
    swap = IborSwap(
        settlement_date,
        "12Y",
        swap_type,
        6.10 / 100.0,
        fixed_freq_type,
        fixed_dcc_type,
        cal_type=cal,
    )
    swaps.append(swap)
    swap = IborSwap(
        settlement_date,
        "15Y",
        swap_type,
        5.90 / 100.0,
        fixed_freq_type,
        fixed_dcc_type,
        cal_type=cal,
    )
    swaps.append(swap)
    swap = IborSwap(
        settlement_date,
        "20Y",
        swap_type,
        5.60 / 100.0,
        fixed_freq_type,
        fixed_dcc_type,
        cal_type=cal,
    )
    swaps.append(swap)
    swap = IborSwap(
        settlement_date,
        "25Y",
        swap_type,
        5.55 / 100.0,
        fixed_freq_type,
        fixed_dcc_type,
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
        check_refit_flag=False,
        do_build=do_build,
        **optional_interp_params,
    )

    calibrator = IborSingleCurveSmoothingCalibrator(init_curve)

    # here we go
    smooth_param = 1e-2
    curve, report = calibrator.fit(
        smoothness=smooth_param, report_progress=REPORT_PROGRESS
    )
    if REPORT_PROGRESS:
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


########################################################################################

if __name__ == "__main__":
    test_smooth_fit_simple(InterpTypes.LINEAR_ZERO_RATES)
    test_smooth_fit(InterpTypes.FLAT_FWD_RATES)
    # test_smooth_fit(InterpTypes.FINCUBIC_ZERO_RATES)
