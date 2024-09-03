import numpy as np
import pandas as pd
from typing import List, Union


from ...utils.tenor import Tenor
from ...utils.date import Date, datediff
from ...utils.day_count import DayCountTypes
from ...utils.global_vars import g_basis_point
from ...market.curves.discount_curve import DiscountCurve
from ...market.curves.discount_curve_pwf_onf import DiscountCurvePWFONF
from ...market.curves.composite_discount_curve import CompositeDiscountCurve

from ...products.rates.ibor_single_curve import IborSingleCurve
from ...products.rates.ibor_single_curve_par_shocker import (
    IborSingleCurveParShocker,
)

PV_PREFIX = "PV:"
DV01_PREFIX = "DV01:"
ROLL_PREFIX = "ROLL:"


def par_rate_risk_report(
    base_curve: IborSingleCurve,
    trades: list,
    trade_labels: list = None,
    bump_size=1.0 * g_basis_point,
):
    """Calculate deltas (change in value to 1bp bump) of the trades to all
    benchmarks in the base curve. Supported trades are depos, fras, swaps.
    trade_labels are used to identify trades in the output, if not provided
    simple ones are generated

    Args:
        base_curve (IborSingleCurve): Base curve to be bumped
        trades (list): a list of trades to calculate deltas of
        trade_labels (list, optional): trade labels to identify trades in
        the output. Defaults to None in which case these are auto generted
        bump_size (float, optional): How big of a bump to apply to bechmarks.
        Output always expressed as change in value per 1 bp.
        Defaults to 1.0*g_basis_point.

    Returns:
        (base_values, risk_report):
            base_value is a dictionary with trade_labels as keys and base
            trade values as values risk_report is a dataframe with benchmarks
            as rows and a column of par deltas per trade, and a total for
            all trades
    """
    curve_shocker = IborSingleCurveParShocker(base_curve)
    benchmarks_report = curve_shocker.benchmarks_report()
    n_benchmarks = curve_shocker.n_benchmarks()
    n_trades = len(trades)

    risk_report = benchmarks_report[
        ["type", "start_date", "maturity_date", "market_rate"]
    ].copy()

    if trade_labels is None:
        trade_labels = [f"trade_{n:03d}" for n in range(n_trades)]

    base_values = {}
    for trade, trade_label in zip(trades, trade_labels):
        base_values[trade_label] = trade.value(base_curve.value_dt, base_curve)

    for benchmark_idx in range(n_benchmarks):
        bumped_curve = curve_shocker.apply_bump_to_benchmark(
            benchmark_idx, bump_size
        )

        for trade_idx, trade in enumerate(trades):
            trade_label = trade_labels[trade_idx]
            base_value = base_values[trade_label]
            bumped_value = trade.value(bumped_curve.value_dt, bumped_curve)
            par_delta = (bumped_value - base_value) / bump_size * g_basis_point
            risk_report.loc[benchmark_idx, trade_label] = par_delta

    risk_report["total"] = risk_report[trade_labels].sum(axis=1)
    return base_values, risk_report


def forward_rate_risk_report(
    base_curve: DiscountCurve,
    grid_last_date: Date,
    grid_bucket_tenor: str,
    trades: list,
    trade_labels: list = None,
    bump_size=1.0 * g_basis_point,
):
    """Generate forward rate deltas (forward delta ladder) risk report, which
    is the sensitivity of trades to bucketed shocks of the instantaneous (ON)
    forward rates. Here shock_i is applied to the ON forward rates over
    [t_i, t_{i+1}] where {t_i} is a time grid from base_curve.valuation_date
    to grid_last_date with buckets of size grid_bucket_tenor

    Args:
        base_curve (DiscountCurve): base curve to apply bumps to
        grid_last_date: the last date for the grid that defines shocks
        grid_bucket_tenor (str): spacing of grid points as tenor string eg '3M'
        trades (list): a list of trades to calculate deltas of
        trade_labels (list, optional): trade labels to identify trades in
        output. Defaults to None in which case these are auto generted
        bump_size (float, optional): How big of a bump to apply to bechmarks.
        Output always expressed as change in value per 1 bp.
        Defaults to 1.0*g_basis_point.

    Returns:
        (dict, Dataframe): (base_values, risk_report)
            base_value is a dictionary with trade_labels as keys and base
            trade values as values
            risk_report is a dataframe with bump details for rows and a
            column of forward rate deltas per trade, and a total for all trades

    """
    valuation_date = base_curve.value_dt
    grid, grid_labels = _grid_from_dates_tenor(
        grid_last_date, grid_bucket_tenor, valuation_date
    )

    base_values, risk_report, *_ = forward_rate_risk_report_custom_grid(
        base_curve, grid, trades, grid_labels, trade_labels, bump_size
    )
    return base_values, risk_report


def forward_rate_risk_report_custom_grid(
    base_curve: DiscountCurve,
    grid: List[Date],
    trades: list,
    grid_labels: list = None,
    trade_labels: list = None,
    bump_size=1.0 * g_basis_point,
):
    """Generate forward rate deltas risk report, which is the sensitivity of
    trades to bucketed shocks of the instantaneous (ON) forward rates. Here
    shock_i is applied to the ON forward rates
    over [t_i, t_{i+1}] where {t_i} is the 'grid' argument

    Args:
        base_curve (DiscountCurve): base curve to apply bumps to
        grid: Date grid that defines ON forward rate bumps. Note that
        curve.valuation_date must be explictly included (if so desired)
        trades (list): a list of trades to calculate deltas of
        trade_labels (list, optional): trade labels to identify trades in the
        output. Defaults to None in which case these are auto generted
        bump_size (float, optional): How big of a bump to apply to bechmarks.
        Output always expressed as change in value per 1 bp.
        Defaults to 1.0*g_basis_point.

    Returns:
        (dict, Dataframe): (base_values, risk_report)
            base_value is a dictionary with trade_labels as keys and base
            trade values as values
            risk_report is a dataframe with bump details for rows and a
            column of forward rate deltas per trade, and a total for all trades

    """
    risk_report = pd.DataFrame(columns=["type", "start_date", "maturity_date"])
    risk_report["start_date"] = grid[:-1]
    risk_report["maturity_date"] = grid[1:]
    risk_report["type"] = "IborFRA"

    asof = base_curve.value_dt
    if grid_labels is None:
        start_in_days = [datediff(asof, d) for d in risk_report["start_date"]]
        tenor_in_days = [
            datediff(s, m)
            for s, m in zip(
                risk_report["start_date"], risk_report["maturity_date"]
            )
        ]
        grid_labels = [
            f"{sd}Dx{td}D" for sd, td in zip(start_in_days, tenor_in_days)
        ]
    risk_report["bucket_label"] = grid_labels

    n_trades = len(trades)
    if trade_labels is None:
        trade_labels = [f"trade_{n:03d}" for n in range(n_trades)]

    base_values = {}
    first_period_carry = {}
    for trade, trade_label in zip(trades, trade_labels):
        val_res = trade.value(asof, base_curve, pv_only=False)
        base_values[trade_label] = val_res[0]

        cf_report = val_res[1]
        cf_report_first_period = cf_report[
            (cf_report["payment_date"] > grid[0])
            & (cf_report["payment_date"] <= grid[1])
        ]
        first_period_carry[trade_label] = cf_report_first_period[
            "payment_pv"
        ].sum()

    for fwdrate_idx in range(len(risk_report)):
        start_date = risk_report.loc[fwdrate_idx, "start_date"]
        maturity_date = risk_report.loc[fwdrate_idx, "maturity_date"]
        base_rate = base_curve.fwd_rate(
            start_date, maturity_date, DayCountTypes.SIMPLE
        )
        risk_report.loc[fwdrate_idx, "market_rate"] = base_rate

        fwd_rate_shock = DiscountCurvePWFONF.brick_wall_curve(
            base_curve.value_dt, start_date, maturity_date, bump_size
        )
        bumped_curve = CompositeDiscountCurve([base_curve, fwd_rate_shock])

        for trade_idx, trade in enumerate(trades):
            trade_label = trade_labels[trade_idx]
            base_value = base_values[trade_label]
            bumped_value = trade.value(bumped_curve.value_dt, bumped_curve)
            par_delta = (bumped_value - base_value) / bump_size * g_basis_point
            risk_report.loc[fwdrate_idx, DV01_PREFIX + trade_label] = par_delta

    risk_report[DV01_PREFIX + "total"] = risk_report[
        [DV01_PREFIX + l for l in trade_labels]
    ].sum(axis=1)
    return base_values, risk_report, first_period_carry


def carry_rolldown_report(
    base_curve: DiscountCurve,
    grid_last_date: Date,
    grid_bucket_tenor: Union[str, Tenor],
    trades: list,
    trade_labels: list = None,
    bump_size=1.0 * g_basis_point,
):
    """Generate carry and rolldown risk report based on the sensitivity of
    trades to bucketed shocks of the instantaneous (ON) forward rates. Here
    shock_i is applied to the ON forward rates over [t_i, t_{i+1}] where {t_i}
    is a time grid from base_curve.valuation_date to grid_last_date with
    buckets of size grid_bucket_tenor

    Args:
        base_curve (DiscountCurve): base curve to apply bumps to
        grid_last_date: the last date for the grid that defines shocks
        grid_bucket_tenor (str): spacing of grid points as tenor string e.g'3M'
        trades (list): a list of trades to calculate deltas of
        trade_labels (list, optional): trade labels to identify trades in the
        output. Defaults to None in which case these are auto generted
        bump_size (float, optional): How big of a bump to apply to bechmarks.
        Output always expressed as change in value per 1 bp.
        Defaults to 1.0*g_basis_point.

    Returns:
        (dict, Dataframe): (base_values, risk_report)
            base_value is a dictionary with trade_labels as keys and base
            trade values as values
            risk_report is a dataframe with time buickets for rows and
            carry/rolldown and DV01 columns
            per trade, and a total for all trades for each measure.
            Columns marked 'ROLL' have carry in the
            first row and rolldown in all the others

    """
    valuation_date = base_curve.value_dt
    grid, grid_labels = _grid_from_dates_tenor(
        grid_last_date, grid_bucket_tenor, valuation_date
    )

    base_values, risk_report, first_period_carry = (
        forward_rate_risk_report_custom_grid(
            base_curve, grid, trades, grid_labels, trade_labels, bump_size
        )
    )

    trade_labels = list(base_values.keys())
    for label in trade_labels:
        risk_report.loc[0, ROLL_PREFIX + label] = first_period_carry[label]
        rate_change = -risk_report["market_rate"].diff()
        risk_report.loc[1:, ROLL_PREFIX + label] = (
            risk_report.loc[1:, DV01_PREFIX + label]
            * rate_change[1:]
            / g_basis_point
        )

    risk_report[ROLL_PREFIX + "total"] = risk_report[
        [ROLL_PREFIX + l for l in trade_labels]
    ].sum(axis=1)

    return base_values, risk_report, first_period_carry


def parallel_shift_ladder_report(
    base_curve: DiscountCurve,
    curve_shifts: np.ndarray,
    trades: list,
    trade_labels: list = None,
):

    curve_shift_labels = [f"SHIFT:{s/g_basis_point:.1f}" for s in curve_shifts]
    risk_report = pd.DataFrame(
        columns=["shift_label", "shift_bp"],
        data=zip(
            curve_shift_labels,
            curve_shifts / g_basis_point,
        ),
    )

    n_trades = len(trades)
    if trade_labels is None:
        trade_labels = [f"trade_{n:03d}" for n in range(n_trades)]

    base_values = {}
    for trade, trade_label in zip(trades, trade_labels):
        base_values[trade_label] = trade.value(base_curve.value_dt, base_curve)

    for shift_idx, shift in enumerate(curve_shifts):
        fwd_rate_shock = DiscountCurvePWFONF.flat_curve(
            base_curve.value_dt, shift
        )
        bumped_curve = CompositeDiscountCurve([base_curve, fwd_rate_shock])

        for trade_idx, trade in enumerate(trades):
            trade_label = trade_labels[trade_idx]
            bumped_value = trade.value(bumped_curve.value_dt, bumped_curve)
            risk_report.loc[shift_idx, PV_PREFIX + trade_label] = bumped_value

    risk_report[PV_PREFIX + "total"] = risk_report[
        [PV_PREFIX + t for t in trade_labels]
    ].sum(axis=1)
    return base_values, risk_report


def _grid_from_dates_tenor(
    grid_last_date, grid_bucket_tenor: Union[str, Tenor], valuation_date
):
    bucket_tenor = Tenor.as_tenor(grid_bucket_tenor)

    grid = [valuation_date]
    grid_labels = []

    # TODO: hardcoded stuff, revisit
    spot_days = 2
    d = grid[0].add_weekdays(spot_days)

    # the loop is structured so the grid_last_date is captured by last bucket
    count = 0
    while d < grid_last_date:
        d = d.add_tenor(bucket_tenor)
        grid.append(d)
        grid_labels.append(f"{count *bucket_tenor}x{bucket_tenor}")
        count += 1
    grid[-1] = grid_last_date
    return grid, grid_labels
