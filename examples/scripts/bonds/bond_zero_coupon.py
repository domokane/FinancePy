# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path
_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause
_install_double_click_pause()
import os
import pandas as pd

import add_fp_to_path

from financepy.utils.date import Date
from financepy.products.bonds.bond_zero import BondZero
from financepy.products.bonds.bond import YTMCalcType


PLOT_GRAPHS = False

########################################################################################


def test_bond_zero():

    issue_dt = Date(25, 7, 2022)
    maturity_dt = Date(24, 10, 2022)
    face_amount = 100.0
    issue_price = 99.6410

    bond = BondZero(issue_dt, maturity_dt, issue_price)

    settle_dt = Date(8, 8, 2022)
    clean_price = 99.6504

    ytm = bond.yield_to_maturity(settle_dt, clean_price)

    accrued_interest = bond.accrued_interest(settle_dt, face_amount)

    print("YTM", "accrued")
    print(ytm, accrued_interest)


########################################################################################


def test_bond_zero_ror():

    # Directory where *this script* is located
    here = os.path.dirname(os.path.abspath(__file__))

    # Path to your local data file inside a "data" folder
    data_path = os.path.join(here, "data", "test_cases_bond_zero_ror.csv")

    #    path = ".//data//test_cases_bond_zero_ror.csv"
    df = pd.read_csv(data_path, parse_dates=["buy_date", "sell_date"])

    # A 1-year bond with zero coupon per year. code: 092103011
    bond = BondZero(
        issue_dt=Date(23, 7, 2021),
        maturity_dt=Date(24, 8, 2022),
        issue_price=97.67,
    )

    print(
        "bond_code",
        "buy_date",
        "buy_ytm",
        "buy_price",
        "sell_date",
        "sell_ytm",
        "sell_price",
        "simple_return",
        "irr",
    )

    for row in df.itertuples(index=False):

        buy_dt = Date(row.buy_date.day, row.buy_date.month, row.buy_date.year)
        sell_dt = Date(row.sell_date.day, row.sell_date.month, row.sell_date.year)

        buy_price = bond.dirty_price_from_ytm(buy_dt, row.buy_ytm, YTMCalcType.ZERO)
        sell_price = bond.dirty_price_from_ytm(sell_dt, row.sell_ytm, YTMCalcType.ZERO)

        simple, irr, _ = bond.calc_ror(buy_dt, sell_dt, row.buy_ytm, row.sell_ytm)

        print(
            row.bond_code,
            buy_dt,
            row.buy_ytm,
            buy_price,
            sell_dt,
            row.sell_ytm,
            sell_price,
            simple,
            irr,
        )


########################################################################################

test_bond_zero()
test_bond_zero_ror()
