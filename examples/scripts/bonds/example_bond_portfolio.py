# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import os
import datetime as dt
import pandas as pd


from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.bonds.bond import Bond
from financepy.utils.date import Date, from_datetime

# ============================================================================
# FINANCEPY EXAMPLES - Bond
# ============================================================================



########################################################################################




########################################################################################

# ============================================================================
# 1. BOND PORTFOLIO
# ============================================================================
# What this section demonstrates:
# Solves for the yield that reproduces the observed bond price. This checks the inverse relationship between price and yield.
# Calculates coupon interest earned since the previous coupon date and illustrates the clean/dirty price adjustment.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. BOND PORTFOLIO")
print("=" * 78)

path = os.path.join(os.path.dirname(__file__), "./data/gilt_bond_prices.txt")
bond_dataframe = pd.read_csv(path, sep="\t")
bond_dataframe["mid"] = 0.5 * (bond_dataframe["bid"] + bond_dataframe["ask"])

freq_type = FrequencyTypes.SEMI_ANNUAL
dc_type = DayCountTypes.ACT_ACT_ICMA

settle_dt = Date(19, 9, 2012)

print("DCTYPE", "MATDATE", "CPN", "PRICE", "ACCD", "YTM")

for dc_type in DayCountTypes:

    if dc_type == DayCountTypes.ZERO:
        continue
    for _, bond in bond_dataframe.iterrows():

        date_string = bond["maturity"]
        mat_date_time = dt.datetime.strptime(date_string, "%d-%b-%y")
        maturity_dt = from_datetime(mat_date_time)
        issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)
        coupon = bond["coupon"] / 100.0
        clean_price = bond["mid"]

        bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

        ytm = bond.yield_to_maturity(settle_dt, clean_price)
        accrued_interest = bond.accrued_interest(settle_dt, 100.0)

        print(
            dc_type,
            maturity_dt,
            coupon * 100.0,
            clean_price,
            accrued_interest,
            ytm * 100.0,
        )

