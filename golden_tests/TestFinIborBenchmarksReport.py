import sys

sys.path.append("..")

import pandas as pd
from os.path import dirname, join

from financepy.utils.date import Date
from financepy.utils.global_types import SwapTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.calendar import CalendarTypes
from financepy.market.curves.interpolator import InterpTypes
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.rates.ibor_fra import IborFRA
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_single_curve import IborSingleCurve

from financepy.products.rates.ibor_benchmarks_report import (
    ibor_benchmarks_report,
    dataframe_to_benchmarks,
)

from FinTestCases import FinTestCases, globalTestCaseMode

test_cases = FinTestCases(__file__, globalTestCaseMode)


def test_ibor_benchmarks_report():
    valuation_date = Date(6, 10, 2001)
    cal = CalendarTypes.UNITED_KINGDOM
    interp_type = InterpTypes.FLAT_FWD_RATES

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

    # Create but do not build the initial curve
    do_build = True
    curve = IborSingleCurve(
        valuation_date,
        depos,
        fras,
        swaps,
        interp_type,
        check_refit=False,
        do_build=do_build,
    )

    benchmarks_report = ibor_benchmarks_report(curve)

    # print(benchmarks_report)

    # Confirm that there are no NaNs. In particular this means that different types of benchmarks
    # return exactly the same keys, just like we want it, with a couple of exceptions
    assert (
        benchmarks_report.drop(columns=["fixed_freq_type", "fixed_leg_type"])
        .isnull()
        .values.any()
    ) == False


def test_dataframe_to_benchmarks():
    path = dirname(__file__)
    filename = "ibor_benchmarks_example.csv"
    full_filename_path = join(path, "data", filename)

    asof = Date(6, 10, 2001)

    df = pd.read_csv(full_filename_path, index_col=0)

    df["start_date"] = pd.to_datetime(
        df["start_date"], format="%d-%b-%y"
    )  # allow tenors
    df["maturity_date"] = pd.to_datetime(
        df["maturity_date"], format="%d-%b-%y"
    )  # allow tenors

    benchmarks = dataframe_to_benchmarks(
        df, asof_date=asof, calendar_type=CalendarTypes.UNITED_KINGDOM
    )

    assert len(benchmarks["IborDeposit"]) == 2
    assert len(benchmarks["IborFRA"]) == 1
    assert len(benchmarks["IborSwap"]) == 10


try:
    test_ibor_benchmarks_report()
    test_dataframe_to_benchmarks()
except Exception as e:
    print(f"Unexpected error:{e}", sys.exc_info()[0])
