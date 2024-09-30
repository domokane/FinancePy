import sys

sys.path.append("..")

from financepy.utils.global_types import SwapTypes
from financepy.utils.global_vars import g_basis_point
from financepy.market.curves.interpolator import InterpTypes
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.rates.ibor_fra import IborFRA
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_single_curve_par_shocker import (
    IborSingleCurveParShocker,
)
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
from financepy.utils.calendar import CalendarTypes

from FinTestCases import FinTestCases, globalTestCaseMode

test_cases = FinTestCases(__file__, globalTestCaseMode)


def test_ibor_curve_par_rate_shocker():
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

    base_curve = IborSingleCurve(
        valuation_date,
        depos,
        fras,
        swaps,
        InterpTypes.FLAT_FWD_RATES,
    )
    curve_shocker = IborSingleCurveParShocker(base_curve)
    mat_dates = curve_shocker.benchmarks_report()["maturity_date"].values

    # size of bump
    par_rate_bump = 1 * g_basis_point

    # expected forward rate changes in the periods before and after the maturity date of the bumped benchmark
    # in basis points
    expected_fwd_rate_changes = {
        1: (1.00000, 0.0),
        2: (1.00000, -0.50715),
        4: (2.06214, -2.16773),
    }

    # Which benchmarks we test-bump. They are
    # - The first (and only) depo we used in construction. Here index = 1 (not 0) because
    # ibor_single_curve creates a synthetic deposit for the stub implictly
    # - The only FRA
    # - the 2Y swap
    benchmark_idxs = [1, 2, 4]
    for benchmark_idx in benchmark_idxs:
        bumped_curve = curve_shocker.apply_bump_to_benchmark(
            benchmark_idx, par_rate_bump
        )

        d1 = mat_dates[benchmark_idx - 1]
        d2 = mat_dates[benchmark_idx]
        d3 = mat_dates[benchmark_idx + 1]

        base_fwd_before = base_curve.fwd_rate(d1, d2)
        base_fwd_after = base_curve.fwd_rate(d2, d3)
        bumped_fwd_before = bumped_curve.fwd_rate(d1, d2)
        bumped_fwd_after = bumped_curve.fwd_rate(d2, d3)

        actual_fwd_rate_changes = (
            (bumped_fwd_before - base_fwd_before) / g_basis_point,
            (bumped_fwd_after - base_fwd_after) / g_basis_point,
        )

        assert round(actual_fwd_rate_changes[0], 3) == round(
            expected_fwd_rate_changes[benchmark_idx][0], 3
        )
        assert round(actual_fwd_rate_changes[1], 3) == round(
            expected_fwd_rate_changes[benchmark_idx][1], 3
        )


test_ibor_curve_par_rate_shocker()
