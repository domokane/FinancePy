# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import sys, os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from financepy.products.bonds.bond_zero_curve import BondZeroCurve
from financepy.products.bonds.bond import Bond
from financepy.utils.date import Date, from_datetime
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
import datetime as dt
import os


import pandas as pd

path = os.path.join(os.path.dirname(__file__), "./data/gilt_bond_prices.txt")
bond_dataframe = pd.read_csv(path, sep="\t")
bond_dataframe["mid"] = 0.5 * (bond_dataframe["bid"] + bond_dataframe["ask"])

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

bond_curve = BondZeroCurve(settle_dt, bonds, clean_prices)

########################################################################################


def test_zero_curve():

    maturity_dt = Date(7, 3, 2013)
    zero_rate = bond_curve.zero_rate(maturity_dt)
    assert round(zero_rate, 4) == 0.0022

    maturity_dt = Date(7, 9, 2019)
    zero_rate = bond_curve.zero_rate(maturity_dt)
    assert round(zero_rate, 4) == 0.0128

    maturity_dt = Date(22, 1, 2060)
    zero_rate = bond_curve.zero_rate(maturity_dt)
    assert round(zero_rate, 4) == 0.0351
