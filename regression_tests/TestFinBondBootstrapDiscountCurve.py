# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import os
import datetime as dt
import pandas as pd
import time
import add_fp_to_path

from financepy.products.bonds.bond import Bond
from financepy.products.bonds import get_bond_market_conventions, BondMarkets
from financepy.market.curves import BondBootstrapDiscountCurve
from financepy.market.curves.interpolator import InterpTypes
from financepy.utils.date import Date, from_datetime
import matplotlib.pyplot as plt

from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

PLOT_GRAPHS = False

########################################################################################

def test_bond_bootstrap_discount_curve():

    path = os.path.join(os.path.dirname(__file__), "./data/gilt_bond_prices.txt")
    bond_dataframe = pd.read_csv(path, sep="\t")
    bond_dataframe["mid"] = 0.5 * (bond_dataframe["bid"] + bond_dataframe["ask"])

    acc_dc_type, freq_type, spot_days, ex_div_days, cal = get_bond_market_conventions(
        BondMarkets.UNITED_KINGDOM
    )

    today = Date(18, 9, 2012)
    settle_dt = today.add_weekdays(spot_days)

    bonds = []
    clean_prices = []

    for _, bond_row in bond_dataframe.iterrows():
        date_string = bond_row["maturity"]
        mat_date_time = dt.datetime.strptime(date_string, "%d-%b-%y")
        maturity_dt = from_datetime(mat_date_time)
        issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)
        coupon = bond_row["coupon"] / 100.0
        clean_price = bond_row["mid"]

        bond = Bond(
            issue_dt,
            maturity_dt,
            coupon,
            freq_type,
            acc_dc_type,
            ex_div_days=ex_div_days,
            cal_type=cal,
        )

        bonds.append(bond)
        clean_prices.append(clean_price)

    for interp_type in InterpTypes:

        start = time.perf_counter()

        bond_curve = BondBootstrapDiscountCurve(
            settle_dt, bonds, clean_prices, interp_type, check_refit_flag=True
        )

        elapsed = time.perf_counter() - start

        test_cases.header("INTERPOLATION", "DATE", "ZERO RATE",
                          "DISCOUNT_FACTOR", "ELAPSED TIME")

        for _, bond_row in bond_dataframe.iterrows():

            date_string = bond_row["maturity"]
            mat_date_time = dt.datetime.strptime(date_string, "%d-%b-%y")
            maturity_dt = from_datetime(mat_date_time)
            zero_rate = bond_curve.zero_rate(maturity_dt)
            df = bond_curve.df(maturity_dt)
            test_cases.print(interp_type, maturity_dt, zero_rate*100, df, elapsed)

        if PLOT_GRAPHS:
            bond_curve.plot_zero_rates("BOND EXACT CURVE - " + str(interp_type))
            bond_curve.plot_fwd_rates("BOND EXACT CURVE - " + str(interp_type))
            plt.show()

########################################################################################

test_bond_bootstrap_discount_curve()

test_cases.compare_test_cases()
