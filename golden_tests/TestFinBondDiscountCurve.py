# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import os

import datetime as dt
import pandas as pd

import add_fp_to_path

from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date, from_datetime

from financepy.products.bonds.bond import Bond
from financepy.products.bonds.bond_exact_zero_curve import BondExactZeroCurve


from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_bond_discount_curve():

    path = os.path.join(os.path.dirname(__file__),
                        ".//data//gilt_bond_prices.txt")
    bond_dataframe = pd.read_csv(path, sep="\t")
    bond_dataframe["mid"] = 0.5 * \
        (bond_dataframe["bid"] + bond_dataframe["ask"])

    freq_type = FrequencyTypes.SEMI_ANNUAL
    settle_dt = Date(19, 9, 2012)
    ex_div_days = 0

    dc_type = DayCountTypes.ACT_ACT_ICMA

    bonds = []
    clean_prices = []
    dirty_prices = []

    for _, bond in bond_dataframe.iterrows():
        date_string = bond["maturity"]
        mat_dt_time = dt.datetime.strptime(date_string, "%d-%b-%y")
        maturity_dt = from_datetime(mat_dt_time)
        issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)

        coupon = bond["coupon"] / 100.0
        clean_price = bond["mid"]

        clean_prices.append(clean_price)

        bond = Bond(issue_dt, maturity_dt, coupon,
                    freq_type, dc_type, ex_div_days)

        accrued = bond.accrued_interest(settle_dt)

        dirty_price = clean_price + accrued

        bonds.append(bond)
        dirty_prices.append(dirty_price)

        ytm = bond.yield_to_maturity(settle_dt, clean_price)
        accrued_int = bond.accrued_int
        accd_days = bond.accrued_days

        bond_curve = BondExactZeroCurve(settle_dt,
                                        bonds,
                                        clean_prices
                                        )

        ###############################################################################


test_bond_discount_curve()

test_cases.compare_test_cases()
