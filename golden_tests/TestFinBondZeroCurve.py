###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import os
import sys
import datetime as dt

sys.path.append("..")

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.bonds.bond import Bond
from financepy.products.bonds.bond_zero_curve import BondZeroCurve
from financepy.utils.date import Date, from_datetime
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes

test_cases = FinTestCases(__file__, globalTestCaseMode)

plotGraphs = False

###############################################################################


def test_BondZeroCurve():

    import pandas as pd

    path = os.path.join(
        os.path.dirname(__file__), "./data/gilt_bond_prices.txt"
    )
    bond_dataframe = pd.read_csv(path, sep="\t")
    bond_dataframe["mid"] = 0.5 * (
        bond_dataframe["bid"] + bond_dataframe["ask"]
    )

    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    settlement = Date(19, 9, 2012)

    bonds = []
    clean_prices = []

    for _, bondRow in bond_dataframe.iterrows():
        date_string = bondRow["maturity"]
        mat_date_time = dt.datetime.strptime(date_string, "%d-%b-%y")
        maturity_dt = from_datetime(mat_date_time)
        issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)
        coupon = bondRow["coupon"] / 100.0
        clean_price = bondRow["mid"]
        bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)
        bonds.append(bond)
        clean_prices.append(clean_price)

    ###########################################################################

    bondCurve = BondZeroCurve(settlement, bonds, clean_prices)

    test_cases.header("DATE", "ZERO RATE")

    for _, bond in bond_dataframe.iterrows():

        date_string = bond["maturity"]
        mat_date_time = dt.datetime.strptime(date_string, "%d-%b-%y")
        maturity_dt = from_datetime(mat_date_time)
        zero_rate = bondCurve.zero_rate(maturity_dt)
        test_cases.print(maturity_dt, zero_rate)

    if plotGraphs:
        bondCurve.plot("BOND CURVE")


###############################################################################

test_BondZeroCurve()
test_cases.compareTestCases()
