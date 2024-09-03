import pandas as pd
from typing import Union
from datetime import datetime


from ...utils.error import FinError
from ...utils.global_types import SwapTypes
from ...utils.calendar import CalendarTypes
from ...utils.day_count import DayCountTypes
from ...utils.date import Date, from_datetime
from ...utils.frequency import FrequencyTypes
from ...market.curves.discount_curve import DiscountCurve
from ...products.rates.ibor_single_curve import IborSingleCurve

from ...products.rates.ibor_fra import IborFRA
from ...products.rates.ibor_swap import IborSwap
from ...products.rates.ibor_deposit import IborDeposit

def benchmarks_report(
    benchmarks,
    valuation_date: Date,
    discount_curve: DiscountCurve,
    index_curve: DiscountCurve = None,
    include_objects=False,
):
    """
    Generate a DataFrame with one row per bechmark. A benchmark is any object that has a function
    valuation_details(...) that returns a dictionary of the right shape. Allowed benchmarks at the moment
    are depos, fras and swaps. Various useful information is reported. This is a bit slow
    so do not use in performance-critical spots
    """

    # benchmarks = depos + fras + swaps
    df_bmi = None
    for benchmark in benchmarks:
        res = benchmark.valuation_details(
            valuation_date, discount_curve, index_curve
        )

        if df_bmi is None:
            df_bmi = pd.DataFrame.from_dict(res, orient="index").T
        else:
            df_bmi = df_bmi._append(res, ignore_index=True)

    if include_objects:
        df_bmi["benchmark_objects"] = benchmarks
    return df_bmi


def ibor_benchmarks_report(iborCurve: IborSingleCurve, include_objects=False):
    """
    Generate a DataFrame with one row per bechmark used in constructing a given iborCurve.
    Various useful information is reported. This is a bit slow so do not use in performance-critical spots
    """

    benchmarks = (
        iborCurve.used_deposits + iborCurve.used_fras + iborCurve.used_swaps
    )
    return benchmarks_report(
        benchmarks,
        iborCurve.value_dt,
        iborCurve,
        include_objects=include_objects,
    )


def _date_or_tenor_to_date(
    date_or_tenor: Union[Date, str, datetime], asof_date: Date
):
    if isinstance(date_or_tenor, str):
        return asof_date.add_tenor(date_or_tenor)
    if isinstance(date_or_tenor, Date):
        return date_or_tenor
    if isinstance(date_or_tenor, datetime):
        return from_datetime(date_or_tenor)

    raise FinError(
        f"{date_or_tenor} is of type {type(date_or_tenor)}, expecting a Date or a tenor string or a datetime"
    )


def _date_or_tenor_to_date_or_tenor(date_or_tenor: Union[Date, str, datetime]):
    if isinstance(date_or_tenor, str):
        return date_or_tenor
    if isinstance(date_or_tenor, Date):
        return date_or_tenor
    if isinstance(date_or_tenor, datetime):
        return from_datetime(date_or_tenor)

    raise FinError(
        f"{date_or_tenor} is of type {type(date_or_tenor)}, expecting a Date or a tenor string or a datetime"
    )


def _deposit_from_df_row(
    row: pd.Series, asof_date: Date, calendar_type: CalendarTypes
):
    cls = globals()[row["type"]]
    return cls(
        start_dt=_date_or_tenor_to_date(row["start_date"], asof_date),
        maturity_dt_or_tenor=_date_or_tenor_to_date_or_tenor(
            row["maturity_date"]
        ),
        deposit_rate=row["contract_rate"],
        dc_type=DayCountTypes[row["day_count_type"]],
        notional=row["notional"],
        cal_type=calendar_type,
    )


def _fra_from_df_row(
    row: pd.Series, asof_date: Date, calendar_type: CalendarTypes
):
    cls = globals()[row["type"]]
    return cls(
        start_dt=_date_or_tenor_to_date(row["start_date"], asof_date),
        maturity_dt_or_tenor=_date_or_tenor_to_date_or_tenor(
            row["maturity_date"]
        ),
        fra_rate=row["contract_rate"],
        dc_type=DayCountTypes[row["day_count_type"]],
        notional=row["notional"],
        pay_fixed_rate=SwapTypes[row["fixed_leg_type"]] == SwapTypes.PAY,
        cal_type=calendar_type,
    )


def _swap_from_df_row(
    row: pd.Series, asof_date: Date, calendar_type: CalendarTypes
):
    cls = globals()[row["type"]]
    return cls(
        effective_dt=_date_or_tenor_to_date(row["start_date"], asof_date),
        term_dt_or_tenor=_date_or_tenor_to_date_or_tenor(row["maturity_date"]),
        fixed_leg_type=SwapTypes[row["fixed_leg_type"]],
        fixed_cpn=row["contract_rate"],
        fixed_freq_type=FrequencyTypes[row["fixed_freq_type"]],
        fixed_dc_type=DayCountTypes[row["day_count_type"]],
        notional=row["notional"],
        cal_type=calendar_type,
    )


def _unknown_from_df_row(row: pd.Series):
    instr_type = row["type"]
    raise FinError(f"No benchmark creator found for type {instr_type}")


def dataframe_to_benchmarks(
    df: pd.DataFrame, asof_date: Date, calendar_type: CalendarTypes
):
    """Crete IborBenchmarks from a dataframe. The dataframe should have at least these columns
    with these sample inputs:
                type   start_date maturity_date     day_count_type notional contract_rate fixed_leg_type fixed_freq_type
        0   IborDeposit  06-OCT-2001   09-OCT-2001            ACT_360    100.0         0.042            NaN             NaN
        1       IborFRA  09-JAN-2002   09-APR-2002            ACT_360    100.0         0.042            PAY             NaN
        2      IborSwap  09-OCT-2001   09-OCT-2002  THIRTY_E_360_ISDA  1000000         0.042            PAY     SEMI_ANNUAL

    start_date and maturity_date could be Date, string convertible to Tenor, or datetime
    Args:
        df (pd.DataFrame): dataframe as bove, with benchmark info
        asof_date (Date): if start_date is a tenor string, asof_date is used as an anchor for that
        calendar_type (CalendarTypes): What calendar to use

    Returns:
        dict: Keys are benchmark types as in df['type'], values are lists of benchmarks of that type
    """
    benchmark_creators = {
        "IborDeposit": _deposit_from_df_row,
        "IborFRA": _fra_from_df_row,
        "IborSwap": _swap_from_df_row,
    }
    benchmarks = {}
    for _, row in df.iterrows():
        bm_type = row["type"]
        bm = benchmark_creators.get(bm_type, _unknown_from_df_row)(
            row, asof_date, calendar_type
        )
        benchmarks[bm_type] = benchmarks.get(bm_type, []) + [bm]
    return benchmarks
