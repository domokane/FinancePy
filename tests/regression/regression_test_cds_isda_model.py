#############################################################################
# Copyright (C) 2018-2026 Dominic O'Kane
#############################################################################

"""Official ISDA USD SOFR grid regression (26 April 2021).

Run as a standard FinancePy regression test. The two CSV fixtures belong
in tests/regression/data. FinTestCases records one valuation row per CDS;
the official reference grid is independently checked at full precision.
"""


from FinTestCases import FinTestCases, global_test_case_mode

import csv
from datetime import date, datetime, timedelta
from pathlib import Path
import numpy as np
import pandas as pd

import add_fp_to_path

from financepy.utils.date import Date
from financepy.utils.calendar import CalendarTypes
from financepy.utils.global_types import CDSModelTypes
from financepy.utils.global_vars import CLEAN
from financepy.market.curves.cds_curve import CDSCurve
from financepy.products.credit.cds import CDS
from financepy.utils.calendar import Calendar, BusDayAdjustTypes, DateGenRuleTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.global_types import SwapTypes
from financepy.products.rates.ois import OIS
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.market.curves.ois_curve import OISCurve
from financepy.market.curves.interpolator import InterpTypes
HERE = Path(__file__).resolve().parent
TRADE = date(2021, 4, 26)
N = 10_000_000.0
METRICS = ('upfront', 'cash_settlement', 'clean_price', 'accrued_days', 'accrued_premium')
HEADER = ['Trade Date', 'Maturity Date', 'Currency', 'Coupon', 'Quoted Spread',
          'Recovery', 'Clean Upfront', 'Cash Settlement Amount', 'Clean Price',
          'Days Accrued', 'Accrued Premium', 'Cash Settlement Date', 'Start Date']


verbose = False
FAST = True  # False checks every public pricing wrapper on every grid row.
MONEY_TOL = 0.01
GRID = HERE / 'data' / 'CDS_ISDA_TEST_USD_SOFR_20210426.csv'
REPORT = HERE / 'compare' / 'isda_ois_grid_errors.csv'

for folder in ("golden", "compare", "differences"):
    (HERE / folder).mkdir(exist_ok=True)


test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def as_date(x):
    if isinstance(x, datetime):
        return x.date()
    if isinstance(x, date):
        return x
    if isinstance(x, str) and "-" in x:
        return date.fromisoformat(x)
    n = int(x)
    if n < 100000:
        return date(1899, 12, 30) + timedelta(days=n)
    return datetime.strptime(str(n), '%Y%m%d').date()


def load_inputs(path):
    """Read the original ISDA data exported to two lossless CSV fixtures."""

    path = Path(path)
    rates_path = path.with_name(path.stem + '_rates.csv')
    rates = pd.read_csv(rates_path, float_precision='round_trip')
    grid = pd.read_csv(path, float_precision='round_trip')

    if list(rates.columns) != ['Type', 'Tenor', 'Rates']:
        raise ValueError('Unexpected RFR Curve headers')

    quotes = []
    for kind, tenor, rate in rates.itertuples(index=False, name=None):
        if str(kind).strip() != 'S' or not np.isfinite(rate):
            raise ValueError(f'Invalid swap quote: {kind}, {tenor}, {rate}')
        quotes.append((str(tenor).strip(), float(rate)))

    if len(quotes) != 19 or quotes[0][0] != '1M':
        raise ValueError('Expected 19 USD SOFR swaps starting at 1M')

    if list(grid.columns) != HEADER:
        raise ValueError('Unexpected Test Grid headers')

    cases = []
    for line, row in enumerate(grid.itertuples(index=False, name=None), 2):
        trade, mat, ccy, coupon, spread, rec, up, csa, px, days, ai, settle, step = row

        if as_date(trade) != TRADE or str(ccy).strip() != 'USD':
            raise ValueError(f'Unexpected trade date or currency at row {line}')

        expected = np.array([up, csa, px, days, ai], dtype=float)
        if not np.all(np.isfinite(expected)):
            raise ValueError(f'Missing/nonfinite reference values at row {line}')

        cases.append((line, as_date(mat), float(coupon)/10000., float(spread)/10000.,
                      float(rec), as_date(settle), as_date(step), expected))

    if len(cases) != 7164:
        raise ValueError(f'Expected 7164 grid rows, got {len(cases)}')

    return quotes, cases


def build_ois_curve(quotes):

    value_dt = Date(26, 4, 2021)
    cal = CalendarTypes.WEEKEND
    bd = BusDayAdjustTypes.MODIFIED_FOLLOWING
    spot_dt = Calendar(cal).add_business_days(value_dt, 2)
    dc = DayCountTypes.ACT_360
    swaps = []

    for tenor, rate in quotes:

        swaps.append(OIS(spot_dt, tenor, SwapTypes.PAY, rate,
                         FrequencyTypes.ANNUAL, dc, payment_lag=0,
                         float_freq_type=FrequencyTypes.ANNUAL,
                         float_dc_type=dc, cal_type=cal, bd_type=bd,
                         dg_type=DateGenRuleTypes.BACKWARD))

    # No market deposit in the grid. The first OIS has a single payment,
    # hence Z(spot,T)=1/(1+S*ACT360(spot,T)). Extend its implied flat
    # forward rate back to trade date; encode that DF as a synthetic deposit.
    first_end = swaps[0].maturity_dt
    days = first_end - spot_dt
    spot_days = spot_dt - value_dt
    first_log_df = -np.log1p(quotes[0][1] * days / 360.)

    bridge_rate = np.expm1(-first_log_df * spot_days / days) / (spot_days / 360.)

    deposit = IborDeposit(value_dt, spot_dt, float(bridge_rate), dc,
                          cal_type=cal, bd_type=bd)

    curve = OISCurve(value_dt, [deposit], [], swaps,
                     interp_type=InterpTypes.FLAT_FWD_RATES,
                     check_refit_flag=False, curve_dc_type=DayCountTypes.ACT_365F)

    if verbose:
        # Check directly rather than relying on version-dependent check_refit().
        print(f'Trade={value_dt}, spot={spot_dt}, synthetic bridge rate={bridge_rate:.12g}')
        print(f'{"Tenor":>6} {"Quote (bp)":>14} {"DF":>18} {"Swap PV/N":>16}')

    for (tenor, rate), swap in zip(quotes, swaps):

        residual = swap.value(value_dt, curve) / swap.fixed_leg.notional

        if verbose:
            print(f'{tenor:>6} {rate*10000:14.6f} {curve.df(swap.maturity_dt):18.12f} {residual:16.3e}')

        if not np.isfinite(residual) or abs(residual) > 1e-8:
            raise AssertionError(f'OIS calibration failed: {tenor}, PV/N={residual}')

    if not np.isclose(curve.df(value_dt), 1., atol=1e-13, rtol=0.):
        raise AssertionError('OIS curve is not anchored at trade date')

    return value_dt, curve


def price_grid(value_dt, discount, cases, model_name, report, money_tol,
               fast=False):

    def fd(d):
        return Date(d.day, d.month, d.year)

    model = CDSModelTypes[model_name]

    tol = np.array([money_tol/N, money_tol, 100.*money_tol/N, 1e-8, money_tol])

    # Each conventional spread is a separate constant-hazard calibration.
    # Reuse the curve for both contractual coupons, never across recoveries.

    cache = {}
    failures = {}
    trade_cache = {}
    df_cache = {}
    api_checked = set()

    worst = np.full(len(METRICS), -np.inf)
    worst_row = np.zeros(len(METRICS), dtype=int)
    passed = errors = 0
    regression_rows = []

    with report.open('w', newline='') as f:

        writer = csv.writer(f)

        header = ['excel_row', 'maturity', 'coupon_bp', 'spread_bp', 'recovery',
                  'step_in', 'settlement', 'calibration_clean_pv', 'pass']

        for name in METRICS:

            header += [name+'_isda', name+'_financepy', name+'_error']

        writer.writerow(header + ['exception'])

        for i, (line, mat, coupon, spread, recovery, settle, step, expected) in enumerate(cases, 1):

            key = (mat, spread, recovery, settle, step)
            prefix = [line, mat, coupon*10000, spread*10000, recovery, step, settle]
            residual = np.nan
            actual = np.full(len(METRICS), np.nan)
            ok = False
            exception = ''

            try:
                mat_dt, step_dt, settle_dt = fd(mat), fd(step), fd(settle)

                if key in failures:
                    raise RuntimeError(failures[key])
                if key not in cache:
                    try:
                        benchmark = CDS(step_dt,
                                        mat_dt,
                                        spread,
                                        N,
                                        True,
                                        cal_type=CalendarTypes.WEEKEND)

                        if benchmark.get_settle_dt(value_dt) != settle_dt:
                            raise ValueError('Calibration settlement does not match grid')

                        curve = CDSCurve(value_dt, [benchmark], discount, recovery,
                                         cds_model=model)

                        if curve.cds_model != model:
                            raise AssertionError('Calibration model mismatch')

                        residual = benchmark.value(value_dt, curve, recovery,
                                                   settle_dt=settle_dt)[CLEAN]
                        if not np.isfinite(residual) or abs(residual) > 0.01:
                            raise AssertionError(f'Calibration clean PV={residual}')
                        cache[key] = (curve, residual)

                    except Exception as exc:
                        failures[key] = str(exc)
                        raise

                curve, residual = cache[key]

                # Contractual schedules are independent of spread and recovery.
                trade_key = (mat, step, coupon)
                if trade_key not in trade_cache:
                    trade_cache[trade_key] = CDS(
                        step_dt, mat_dt, coupon, N, True,
                        cal_type=CalendarTypes.WEEKEND,
                    )
                cds = trade_cache[trade_key]

                if settle not in df_cache:
                    df_cache[settle] = discount.df(settle_dt)
                df = df_cache[settle]

                dirty, clean = cds.value(
                    value_dt, curve, recovery, settle_dt=settle_dt
                )

                actual = np.array([
                    clean / (N * df),
                    dirty / df,
                    100.0 * (1.0 - clean / (N * df)),
                    cds.accrued_days(value_dt),
                    cds.accrued_interest(value_dt),
                ], dtype=float)

                # Full mode checks every public method. Fast mode checks
                # the first instance of each maturity/coupon/recovery group
                # and the last case, in addition to every grid comparison.
                api_key = (mat, coupon, recovery, step, settle)
                if not fast or api_key not in api_checked or i == len(cases):
                    public_values = np.array([
                        cds.upfront(value_dt, curve, recovery, settle_dt=settle_dt),
                        cds.cash_settlement_amount(value_dt, curve, recovery, settle_dt=settle_dt),
                        cds.clean_price(value_dt, curve, recovery, settle_dt=settle_dt),
                    ], dtype=float)
                    for public, derived, atol in zip(
                        public_values, actual[:3], (1e-12, 1e-6, 1e-10)
                    ):
                        np.testing.assert_allclose(
                            public, derived, rtol=0.0, atol=atol
                        )
                    api_checked.add(api_key)
                    actual[:3] = public_values

                if not np.all(np.isfinite(actual)):
                    raise ValueError('Nonfinite pricing output')

                # Internal identities are independent of the official grid.
                np.testing.assert_allclose([clean-dirty, actual[1], actual[2]],
                                           [actual[4]*df, N*actual[0]-actual[4], 100.*(1.-actual[0])],
                                           rtol=0., atol=1e-6)
                error = actual - expected

                metric_ok = np.abs(error) <= tol
                ok = bool(np.all(metric_ok))

                for j in range(len(METRICS)):
                    if abs(error[j]) > worst[j]:
                        worst[j] = abs(error[j])
                        worst_row[j] = line

            except Exception as exc:
                errors += 1
                exception = f'{type(exc).__name__}: {exc}'

            passed += int(ok)

            # One regression output row per official CDS case, in grid order.
            # Convert NumPy scalars to Python floats for FinTestCases.
            regression_rows.append((
                int(line), str(mat), float(coupon * 10000),
                float(spread * 10000), float(recovery),
                str(step), str(settle),
                float(actual[0] * N), float(actual[1]),
                float(actual[2]), float(actual[3]), float(actual[4]),
            ))

            output = prefix + [residual, ok]

            for e, a in zip(expected, actual):
                output += [e, a, a-e]

            writer.writerow(output + [exception])
            if verbose and (i % 500 == 0 or i == len(cases)):
                print(f'{i}/{len(cases)} priced; {passed} passed, {errors} exceptions', flush=True)

    if verbose:
        print(f'\nModel: {model_name}; calibrated curves: {len(cache)}')
        print(f'Passed {passed}/{len(cases)}; exceptions {errors}; report {report}')
        print(f'{"Measure":<20} {"Max abs error":>18} {"Tolerance":>18} {"Excel row":>10}')

        for name, error, tolerance, line in zip(METRICS, worst, tol, worst_row):
            if np.isfinite(error):
                print(f'{name:<20} {error:18.10g} {tolerance:18.10g} {line:10d}')
            else:
                print(f'{name:<20} {"N/A":>18} {tolerance:18.10g} {"N/A":>10}')

    if errors == len(cases):
        print('No successful pricing comparisons; inspect the exception column.')
    elif verbose:
        print(f'Max upfront error in dollars (successful comparisons): {worst[0]*N:.10g}')

    ok = passed == len(cases) and len(cache) == 3582
    return ok, regression_rows


########################################################################################


def test_isda_usd_sofr_grid():
    """Build the OIS curve and validate every official CDS grid row."""

    quotes, cases = load_inputs(GRID)
    value_dt, discount = build_ois_curve(quotes)

    REPORT.parent.mkdir(parents=True, exist_ok=True)

    ok, regression_rows = price_grid(
        value_dt, discount, cases, "ISDA", REPORT, MONEY_TOL, fast=FAST
    )

    # The official grid is the independent reference, not the golden file.
    if not ok:
        raise AssertionError(f"ISDA USD SOFR grid regression failed. See {REPORT}")

    test_cases.header(
        "EXCEL_ROW", "MATURITY", "COUPON_BP", "SPREAD_BP", "RECOVERY",
        "STEP_IN", "SETTLEMENT", "UPFRONT_USD", "CASH_SETTLEMENT_USD",
        "CLEAN_PRICE", "ACCRUED_DAYS", "ACCRUED_PREMIUM_USD",
    )

    for row in regression_rows:
        test_cases.print(*row)


########################################################################################


test_isda_usd_sofr_grid()
test_cases.compare_test_cases()
