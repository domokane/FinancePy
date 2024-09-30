import numpy as np
import matplotlib.pyplot as plt

import sys

sys.path.append("..")

from helpers import buildIborSingleCurve
from financepy.utils.date import Date
from financepy.utils.global_vars import g_basis_point
from financepy.utils.global_types import SwapTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.market.curves.interpolator import InterpTypes
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_fra import IborFRA
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.rates.ibor_single_curve import IborSingleCurve
import financepy.products.rates.ibor_curve_risk_engine as re

from FinTestCases import FinTestCases, globalTestCaseMode

test_cases = FinTestCases(__file__, globalTestCaseMode)

# when set to True this file can be run standalone and will produce some useful output.
# Set to False to use as part of a testing framework
DIAGNOSTICS_MODE = False


def test_par_rate_risk_report_cubic_zero():
    valuation_date = Date(6, 10, 2001)
    cal = CalendarTypes.UNITED_KINGDOM
    interp_type = InterpTypes.FINCUBIC_ZERO_RATES

    depoDCCType = DayCountTypes.ACT_360
    fraDCCType = DayCountTypes.ACT_360
    swapType = SwapTypes.PAY
    fixedDCCType = DayCountTypes.THIRTY_E_360_ISDA
    fixed_freqType = FrequencyTypes.SEMI_ANNUAL

    settlement_date, base_curve = _generate_base_curve(
        valuation_date,
        cal,
        interp_type,
        depoDCCType,
        fraDCCType,
        swapType,
        fixedDCCType,
        fixed_freqType,
    )
    trades = _generate_trades(
        valuation_date,
        cal,
        swapType,
        fixedDCCType,
        fixed_freqType,
        settlement_date,
        base_curve,
    )

    # size of bump to apply. In all cases par risk is reported as change in value to 1 bp rate bump
    par_rate_bump = 1 * g_basis_point

    # run the report
    base_values, risk_report = re.par_rate_risk_report(
        base_curve, trades, bump_size=par_rate_bump
    )

    expected_totals = [
        0.00122854,
        -0.25323828,
        -0.24271177,
        -0.01423219,
        0.31617136,
        4.0262114,
        2.03409619,
        -0.3957559,
    ]
    actual_totals = risk_report["total"].values

    if DIAGNOSTICS_MODE:
        trade_labels = list(base_values.keys())
        np.set_printoptions(suppress=True)
        print(base_values)
        print(risk_report["total"].values)
        print(risk_report[trade_labels + ["total"]].sum(axis=0))

    assert max(np.abs(actual_totals - expected_totals)) <= 1e-4


def test_par_rate_risk_report_flat_forward():
    valuation_date = Date(6, 10, 2022)
    base_curve = buildIborSingleCurve(valuation_date, "10Y")
    settlement_date = base_curve.used_swaps[0].effective_dt
    cal = base_curve.used_swaps[0].fixed_leg.cal_type
    fixed_day_count = base_curve.used_swaps[0].fixed_leg.dc_type
    fixed_freq_type = base_curve.used_swaps[0].fixed_leg.freq_type

    trades = _generate_trades(
        valuation_date,
        cal,
        SwapTypes.PAY,
        fixed_day_count,
        fixed_freq_type,
        settlement_date,
        base_curve,
    )

    # size of bump to apply. In all cases par risk is reported as change in value to 1 bp rate bump
    par_rate_bump = 1 * g_basis_point

    # run the report
    base_values, risk_report = re.par_rate_risk_report(
        base_curve, trades, bump_size=par_rate_bump
    )

    expected_totals = [
        -0.08629015,
        -0.20597528,
        -0.08628776,
        -0.07793533,
        -0.05012633,
        -0.00005542,
        -0.00005726,
        -0.00005542,
        -0.00005726,
        -0.00005726,
        -0.00008393,
        -0.00013093,
        -0.0001267,
        1.01065078,
        1.49626527,
        3.99996586,
        0,
        0,
        0,
        0,
        0,
        0,
    ]
    actual_totals = risk_report["total"].values

    if DIAGNOSTICS_MODE:
        trade_labels = list(base_values.keys())
        np.set_printoptions(suppress=True)
        print(base_values)
        print(risk_report["total"].values)
        print(risk_report[trade_labels + ["total"]].sum(axis=0))

    assert max(np.abs(actual_totals - expected_totals)) <= 1e-4


def test_forward_rate_risk_report():
    valuation_date = Date(6, 10, 2001)
    cal = CalendarTypes.UNITED_KINGDOM
    interp_type = InterpTypes.FLAT_FWD_RATES

    depoDCCType = DayCountTypes.ACT_360
    fraDCCType = DayCountTypes.ACT_360
    swapType = SwapTypes.PAY
    fixedDCCType = DayCountTypes.THIRTY_E_360_ISDA
    fixed_freqType = FrequencyTypes.SEMI_ANNUAL

    settlement_date, base_curve = _generate_base_curve(
        valuation_date,
        cal,
        interp_type,
        depoDCCType,
        fraDCCType,
        swapType,
        fixedDCCType,
        fixed_freqType,
    )
    trades = _generate_trades(
        valuation_date,
        cal,
        swapType,
        fixedDCCType,
        fixed_freqType,
        settlement_date,
        base_curve,
    )

    # the grid on which we generate the risk report
    grid_bucket = "3M"
    grid_last_date = max(t.maturity_dt for t in trades)

    # size of bump to apply. In all cases par risk is reported as change in value to 1 bp rate bump
    forward_rate_bump = 1 * g_basis_point

    # run the report
    base_values, risk_report = re.forward_rate_risk_report(
        base_curve,
        grid_last_date,
        grid_bucket,
        trades,
        bump_size=forward_rate_bump,
    )

    expected_totals = [
        0.25196322,
        0.24648603,
        0.48750647,
        0.49286362,
        0.48160453,
        0.47113499,
        0.46547237,
        0.47058739,
        0.45980829,
        0.45481044,
        0.23117766,
        0.22408547,
        0.2190641,
        0.21423566,
        0.21169689,
        0.21396794,
        0.0,
    ]
    actual_totals = risk_report[re.DV01_PREFIX + "total"].values

    if DIAGNOSTICS_MODE:
        dv01_trade_labels = [re.DV01_PREFIX + l for l in base_values.keys()]
        np.set_printoptions(suppress=True)
        print(base_values)
        print(risk_report)
        print(risk_report[re.DV01_PREFIX + "total"].values)
        print(
            risk_report[dv01_trade_labels + [re.DV01_PREFIX + "total"]].sum(
                axis=0
            )
        )

    assert max(np.abs(actual_totals - expected_totals)) <= 1e-4


def test_forward_rate_custom_grid_risk_report():
    valuation_date = Date(6, 10, 2001)
    cal = CalendarTypes.UNITED_KINGDOM
    interp_type = InterpTypes.FLAT_FWD_RATES

    depoDCCType = DayCountTypes.ACT_360
    fraDCCType = DayCountTypes.ACT_360
    swapType = SwapTypes.PAY
    fixedDCCType = DayCountTypes.THIRTY_E_360_ISDA
    fixed_freqType = FrequencyTypes.SEMI_ANNUAL

    settlement_date, base_curve = _generate_base_curve(
        valuation_date,
        cal,
        interp_type,
        depoDCCType,
        fraDCCType,
        swapType,
        fixedDCCType,
        fixed_freqType,
    )
    trades = _generate_trades(
        valuation_date,
        cal,
        swapType,
        fixedDCCType,
        fixed_freqType,
        settlement_date,
        base_curve,
    )

    # the grid on which we generate the risk report
    grid = [
        valuation_date,
        valuation_date.add_tenor("3M"),
        valuation_date.add_tenor("15M"),
        valuation_date.add_tenor("10Y"),
    ]

    # size of bump to apply. In all cases par risk is reported as change in value to 1 bp rate bump
    forward_rate_bump = 1 * g_basis_point

    # run the report
    base_values, risk_report, *_ = re.forward_rate_risk_report_custom_grid(
        base_curve, grid, trades, bump_size=forward_rate_bump
    )

    expected_totals = [0.24374713, 1.70089994, 3.65138343]
    actual_totals = risk_report[re.DV01_PREFIX + "total"].values

    if DIAGNOSTICS_MODE:
        dv01_trade_labels = [re.DV01_PREFIX + l for l in base_values.keys()]
        np.set_printoptions(suppress=True)
        print(base_values)
        print(risk_report)
        print(risk_report[re.DV01_PREFIX + "total"].values)
        print(
            risk_report[dv01_trade_labels + [re.DV01_PREFIX + "total"]].sum(
                axis=0
            )
        )

    assert max(np.abs(actual_totals - expected_totals)) <= 1e-4


def test_carry_rolldown_report():
    valuation_date = Date(6, 10, 2001)
    cal = CalendarTypes.UNITED_KINGDOM
    interp_type = InterpTypes.FLAT_FWD_RATES

    depoDCCType = DayCountTypes.ACT_360
    fraDCCType = DayCountTypes.ACT_360
    swapType = SwapTypes.PAY
    fixedDCCType = DayCountTypes.THIRTY_E_360_ISDA
    fixed_freqType = FrequencyTypes.SEMI_ANNUAL

    settlement_date, base_curve = _generate_base_curve(
        valuation_date,
        cal,
        interp_type,
        depoDCCType,
        fraDCCType,
        swapType,
        fixedDCCType,
        fixed_freqType,
    )
    trades = _generate_trades(
        valuation_date,
        cal,
        swapType,
        fixedDCCType,
        fixed_freqType,
        settlement_date,
        base_curve,
    )

    # the grid on which we generate the risk report
    grid_bucket = "6M"
    grid_last_date = max(t.maturity_dt for t in trades)

    # run the report
    base_values, risk_report, *_ = re.carry_rolldown_report(
        base_curve,
        grid_last_date,
        grid_bucket,
        trades,
    )

    if DIAGNOSTICS_MODE:
        roll_trade_labels = [re.ROLL_PREFIX + l for l in base_values.keys()]
        np.set_printoptions(suppress=True)
        print(base_values)
        print(risk_report)
        print(risk_report[re.ROLL_PREFIX + "total"].values)
        print(
            risk_report[roll_trade_labels + [re.ROLL_PREFIX + "total"]].sum(
                axis=0
            )
        )

    expected_totals = [
        -21.07588523,
        16.07402002,
        -27.27637637,
        -0.02469482,
        -104.29563056,
        0.32142506,
        35.98144427,
        0.28310627,
        0.0,
    ]
    actual_totals = risk_report[re.ROLL_PREFIX + "total"].values
    assert max(np.abs(actual_totals - expected_totals)) <= 1e-4


def test_parallel_shift_ladder_report():
    valuation_date = Date(6, 10, 2001)
    cal = CalendarTypes.UNITED_KINGDOM
    interp_type = InterpTypes.FLAT_FWD_RATES

    depoDCCType = DayCountTypes.ACT_360
    fraDCCType = DayCountTypes.ACT_360
    swapType = SwapTypes.PAY
    fixedDCCType = DayCountTypes.THIRTY_E_360_ISDA
    fixed_freqType = FrequencyTypes.SEMI_ANNUAL

    settlement_date, base_curve = _generate_base_curve(
        valuation_date,
        cal,
        interp_type,
        depoDCCType,
        fraDCCType,
        swapType,
        fixedDCCType,
        fixed_freqType,
    )
    trades = _generate_trades(
        valuation_date,
        cal,
        swapType,
        fixedDCCType,
        fixed_freqType,
        settlement_date,
        base_curve,
    )

    # the curve shift grids on which we calculate the PV ladder
    curve_shifts = np.linspace(
        -400 * g_basis_point, 400 * g_basis_point, 17, endpoint=True
    )

    # run the report
    base_values, risk_report = re.parallel_shift_ladder_report(
        base_curve,
        curve_shifts,
        trades,
    )

    if DIAGNOSTICS_MODE:
        pv_trade_labels = [re.PV_PREFIX + l for l in base_values.keys()]
        np.set_printoptions(suppress=True)
        print(base_values)
        print(risk_report)
        print(risk_report[re.PV_PREFIX + "total"].values)
        print(
            risk_report[pv_trade_labels + [re.PV_PREFIX + "total"]].sum(axis=0)
        )

        # risk_report.plot('shift_bp', re.PV_PREFIX + 'total')
        x = risk_report["shift_bp"].values
        y = risk_report[re.PV_PREFIX + "total"].values
        plt.plot(x, y - x * (y[-1] - y[0]) / (x[-1] - x[0]))
        plt.show()

    expected_totals = [
        -2407.30636808,
        -2087.15913624,
        -1772.70446365,
        -1463.84150813,
        -1160.47126898,
        -862.49655251,
        -569.82193811,
        -282.35374508,
        -0.0,
        277.32959522,
        549.7236947,
        817.26933932,
        1080.05198694,
        1338.15554197,
        1591.66238438,
        1840.6533982,
        2085.20799948,
    ]
    actual_totals = risk_report[re.PV_PREFIX + "total"].values
    assert max(np.abs(actual_totals - expected_totals)) <= 1e-4


def _generate_trades(
    valuation_date,
    cal,
    swapType,
    fixedDCCType,
    fixed_freqType,
    settlement_date,
    base_curve,
):
    trade1 = IborSwap(
        settlement_date,
        "4Y",
        swapType,
        4.20 / 100.0,
        fixed_freqType,
        fixedDCCType,
        cal_type=cal,
        notional=10000,
    )
    atm = trade1.swap_rate(valuation_date, base_curve)
    trade1.set_fixed_rate(atm)
    trade2 = IborSwap(
        settlement_date.add_tenor("6M"),
        "2Y",
        swapType,
        4.20 / 100.0,
        fixed_freqType,
        fixedDCCType,
        cal_type=cal,
        notional=10000,
    )
    atm = trade2.swap_rate(valuation_date, base_curve)
    trade2.set_fixed_rate(atm)
    trades = [trade1, trade2]
    return trades


def _generate_base_curve(
    valuation_date,
    cal,
    interp_type,
    depoDCCType,
    fraDCCType,
    swapType,
    fixedDCCType,
    fixed_freqType,
):
    depos = []
    spot_days = 2
    settlement_date = valuation_date.add_weekdays(spot_days)
    depo = IborDeposit(
        settlement_date, "3M", 4.2 / 100.0, depoDCCType, cal_type=cal
    )
    depos.append(depo)

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
        4.70 / 100.0,
        fixed_freqType,
        fixedDCCType,
        cal_type=cal,
    )
    swaps.append(swap)
    swap = IborSwap(
        settlement_date,
        "7Y",
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
        interp_type,
    )

    return settlement_date, base_curve


if 1 == 1:
    test_par_rate_risk_report_cubic_zero()
    test_forward_rate_risk_report()
    test_forward_rate_custom_grid_risk_report()
    test_carry_rolldown_report()
    test_parallel_shift_ladder_report()
