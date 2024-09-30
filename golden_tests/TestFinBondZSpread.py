import os
import pandas as pd
import matplotlib.pyplot as plt

import sys

sys.path.append("..")

from financepy.utils.global_vars import g_percent
from financepy.utils.calendar import CalendarTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date, from_datetime
from financepy.market.curves.interpolator import InterpTypes
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.products.bonds.bond import Bond
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_benchmarks_report import (
    dataframe_to_benchmarks,
)

# Set to True to run this file spandalone and see some useful info
DIAGNOSTICS_MODE = False

from FinTestCases import FinTestCases, globalTestCaseMode

test_cases = FinTestCases(__file__, globalTestCaseMode)


def test_z_spread_flat_curve():

    settlement = Date(19, 9, 2012)
    base_curve = DiscountCurveFlat(settlement, flat_rate=1 * g_percent)
    return _test_z_spread_for_curve(base_curve)


def test_z_spread_actual_curve():

    path = os.path.join(
        os.path.dirname(__file__), "./data/GBP_OIS_20120919.csv"
    )
    dfbm = pd.read_csv(path, index_col=0)

    dfbm["base_date"] = pd.to_datetime(
        dfbm["base_date"], errors="ignore", format="%d/%m/%Y"
    )
    dfbm["start_date"] = pd.to_datetime(
        dfbm["start_date"], errors="ignore", format="%d/%m/%Y"
    )  # allow tenors
    dfbm["maturity_date"] = pd.to_datetime(
        dfbm["maturity_date"], errors="ignore", format="%d/%m/%Y"
    )  # allow tenors

    valuation_date = from_datetime(dfbm.loc[0, "base_date"])
    cal = CalendarTypes.UNITED_KINGDOM
    bms = dataframe_to_benchmarks(
        dfbm, asof_date=valuation_date, calendar_type=cal
    )
    depos = bms["IborDeposit"]
    fras = bms["IborFRA"]
    swaps = bms["IborSwap"]

    fras.sort(key=lambda fra: fra.maturity_dt)
    libor_curve = IborSingleCurve(
        valuation_date, depos, fras, swaps, InterpTypes.LINEAR_ZERO_RATES
    )

    return _test_z_spread_for_curve(libor_curve)


def _test_z_spread_for_curve(base_curve: DiscountCurve):
    path = os.path.join(
        os.path.dirname(__file__), "./data/gilt_bond_prices.txt"
    )
    bond_dataframe = pd.read_csv(path, sep="\t")
    bond_dataframe["mid"] = 0.5 * (
        bond_dataframe["bid"] + bond_dataframe["ask"]
    )

    bond_dataframe["maturity"] = pd.to_datetime(
        bond_dataframe["maturity"], format="%d-%b-%y"
    )
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA

    for bdfIndex, bondRow in bond_dataframe.iterrows():
        matDatetime = bondRow["maturity"]
        maturity_dt = from_datetime(matDatetime)
        issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)
        coupon = bondRow["coupon"] / 100.0
        clean_price = bondRow["mid"]
        bond = Bond(issue_dt, maturity_dt, coupon, freq_type, accrual_type)
        z_spread = bond.z_spread(base_curve.value_dt, clean_price, base_curve)
        asset_swap_spread = bond.asset_swap_spread(
            base_curve.value_dt, clean_price, base_curve
        )
        bond_dataframe.loc[bdfIndex, "z_spread"] = z_spread
        bond_dataframe.loc[bdfIndex, "asset_swap_spread"] = asset_swap_spread

    if DIAGNOSTICS_MODE:
        print(bond_dataframe)
        plt.plot(
            bond_dataframe["maturity"],
            bond_dataframe["gross redemption yield"],
            ".",
            label="yield",
        )
        plt.plot(
            bond_dataframe["maturity"],
            bond_dataframe["z_spread"] * 100,
            ".",
            label="z_spread",
        )
        plt.plot(
            bond_dataframe["maturity"],
            bond_dataframe["asset_swap_spread"] * 100,
            ".",
            label="asset_swap_spread",
        )
        plt.legend(loc="best")
        plt.show()

    assert bond_dataframe["z_spread"].isnull().values.any() == False


if __name__ == "__main__":
    test_z_spread_flat_curve()
    test_z_spread_actual_curve()
