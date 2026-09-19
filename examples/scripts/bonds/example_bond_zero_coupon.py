# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import os
import pandas as pd


from financepy.utils.date import Date
from financepy.products.bonds.bond_zero import BondZero
from financepy.products.bonds.bond import YTMCalcType
import numpy as np
import matplotlib.pyplot as plt

# ============================================================================
# FINANCEPY EXAMPLES - BondZero
# ============================================================================


PLOT_GRAPHS = False

########################################################################################




########################################################################################




########################################################################################

# ============================================================================
# 1. BOND ZERO
# ============================================================================
# What this section demonstrates:
# Solves for the yield that reproduces the observed bond price. This checks the inverse relationship between price and yield.
# Calculates coupon interest earned since the previous coupon date and illustrates the clean/dirty price adjustment.

print("\n" + "=" * 78)
print("1. BOND ZERO")
print("=" * 78)

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

# ============================================================================
# 2. BOND ZERO ROR
# ============================================================================
# What this section demonstrates:
# Prices the bond from a yield including accrued interest. Compare it with clean price to see the effect of accrued coupon.
# Calculates holding-period return, IRR and P&L between purchase and sale dates.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("2. BOND ZERO ROR")
print("=" * 78)

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

# =============================================================================
# 3. VISUALISE ZERO-COUPON BOND PRICE VERSUS YIELD
# =============================================================================
# A zero-coupon bond has only one maturity cash flow, making it a particularly
# clean example of the inverse price/yield relationship.
plot_bond = BondZero(Date(25, 7, 2022), Date(24, 10, 2022), 99.6410)
plot_settle_dt = Date(8, 8, 2022)
plot_yields = np.linspace(0.001, 0.08, 80)
plot_prices = [plot_bond.dirty_price_from_ytm(plot_settle_dt, y, YTMCalcType.ZERO) for y in plot_yields]
plt.figure()
plt.plot(plot_yields * 100.0, plot_prices)
plt.xlabel("Yield to maturity (%)")
plt.ylabel("Dirty price per 100 face")
plt.title("Zero-coupon bond price versus yield")
plt.grid(True)
plt.tight_layout()
plt.show()
