# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import os
import datetime as dt
import pandas as pd

import add_fp_to_path

from financepy.products.bonds.bond import Bond
from financepy.products.bonds.bond_exact_zero_curve import BondExactZeroCurve
from financepy.utils.date import Date, from_datetime
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.market.curves.interpolator import Interpolator, InterpTypes

from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

PLOT_GRAPHS = False

########################################################################################


def test_bond_zero_curve():

    path = os.path.join(
        os.path.dirname(__file__), "./data/gilt_bond_prices.txt"
    )
    bond_dataframe = pd.read_csv(path, sep="\t")
    bond_dataframe["mid"] = 0.5 * (
        bond_dataframe["bid"] + bond_dataframe["ask"]
    )

    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    settle_dt = Date(19, 9, 2012)

    bonds = []
    clean_prices = []

    for _, bond_row in bond_dataframe.iterrows():
        date_string = bond_row["maturity"]
        mat_date_time = dt.datetime.strptime(date_string, "%d-%b-%y")
        maturity_dt = from_datetime(mat_date_time)
        issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)
        coupon = bond_row["coupon"] / 100.0
        clean_price = bond_row["mid"]
        bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)
        bonds.append(bond)
        clean_prices.append(clean_price)

    for interp_type in [
        InterpTypes.FLAT_FWD_RATES,
        InterpTypes.LINEAR_FWD_RATES,
        InterpTypes.LINEAR_ZERO_RATES,
    ]:

        if Interpolator.suitable_for_bootstrap(interp_type):

            bond_curve = BondExactZeroCurve(
                settle_dt, bonds, clean_prices, interp_type
            )

            test_cases.header("DATE", "ZERO RATE")

            for _, bond_row in bond_dataframe.iterrows():

                date_string = bond_row["maturity"]
                mat_date_time = dt.datetime.strptime(date_string, "%d-%b-%y")
                maturity_dt = from_datetime(mat_date_time)
                zero_rate = bond_curve.zero_rate(maturity_dt)
                test_cases.print(maturity_dt, zero_rate)

            if PLOT_GRAPHS:
                bond_curve.plot_zero_rates(
                    "BOND EXACT CURVE - " + str(interp_type)
                )
                bond_curve.plot_fwd_rates(
                    "BOND EXACT CURVE - " + str(interp_type)
                )


########################################################################################

test_bond_zero_curve()
test_cases.compare_test_cases()
