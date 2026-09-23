# ============================================================================
# FINANCEPY EXAMPLES - Bond Futures
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
#
# This module demonstrates the main analytics associated with bond futures
# and their deliverable bond baskets.
#
# The examples cover:
#
#   1. Conversion factors
#   2. Principal and total invoice amounts
#   3. Deliverable-basket analysis
#   4. Yield-to-maturity calculations
#   5. Gross and net basis
#   6. Implied repo rates
#   7. Cheapest-to-deliver (CTD) bond selection
#
# The calculations include examples from published CME material together
# with a comparison using Bloomberg market data.
# ============================================================================

import pandas as pd

from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.format_graphs import set_plot_style

from financepy.products.bonds.bond import Bond, YTMCalcType
from financepy.products.bonds.bond_future import BondFuture

set_plot_style()

# ============================================================================
# 1. CONVERSION FACTORS - MARTELLINI AND PRIAULET
# ============================================================================
#
# Reproduce the bond-futures conversion-factor example from Martellini
# and Priaulet (p. 360).
#
# A conversion factor adjusts for differences between the coupon and
# maturity of a deliverable bond and the notional bond underlying the
# futures contract.
# ============================================================================

freq = FrequencyTypes.SEMI_ANNUAL
basis = DayCountTypes.ACT_ACT_ICMA
issue_dt = Date(15, 2, 2004)

bond1 = Bond(issue_dt, Date(15, 8, 2011), 0.0500, freq, basis)
bond2 = Bond(issue_dt, Date(15, 2, 2011), 0.0500, freq, basis)
bond3 = Bond(issue_dt, Date(15, 8, 2010), 0.0575, freq, basis)
bond4 = Bond(issue_dt, Date(15, 2, 2010), 0.0650, freq, basis)
bond5 = Bond(issue_dt, Date(15, 8, 2009), 0.0600, freq, basis)
bond6 = Bond(issue_dt, Date(15, 5, 2009), 0.0550, freq, basis)
bond7 = Bond(issue_dt, Date(15, 11, 2008), 0.0475, freq, basis)

bonds = [
    bond1,
    bond2,
    bond3,
    bond4,
    bond5,
    bond6,
    bond7,
]

first_delivery_dt = Date(1, 3, 2002)
last_delivery_dt = Date(28, 3, 2002)
contract_size = 100000
contract_cpn = 0.06

bfut = BondFuture(
    "TYH2",
    first_delivery_dt,
    last_delivery_dt,
    contract_size,
    contract_cpn,
)

settle_dt = Date(10, 12, 2001)

# Calculate the exchange conversion factor for each deliverable bond.
print("Bond Maturity", "Coupon", "Conversion Factor")
for bond in bonds:
    cf = bfut.conversion_factor(bond)
    print(bond.maturity_dt, bond.cpn * 100, cf)


# ============================================================================
# 2. CME TREASURY FUTURES INVOICE EXAMPLE
# ============================================================================
#
# Reproduce the Treasury-futures invoice calculations from the CME
# "Understanding Treasury Futures" material:
#
# https://www.cmegroup.com/education/files/understanding-treasury-futures.pdf
#
# The futures price is multiplied by the conversion factor to obtain the
# principal invoice amount. Accrued interest is subsequently included in
# the total invoice amount.
# ============================================================================

freq = FrequencyTypes.SEMI_ANNUAL
basis = DayCountTypes.ACT_ACT_ICMA
issue_dt = Date(15, 2, 2004)

print("EXAMPLE FROM CME")
print("================")

settle_dt = Date(10, 10, 2017)

first_del_dt = Date(1, 12, 2017)
last_del_dt = Date(29, 12, 2017)

fut_size = 100000
fut_coupon = 0.06

bfut = BondFuture(
    "TYZ7",
    first_del_dt,
    last_del_dt,
    fut_size,
    fut_coupon,
)

fut_price = 125.265625

# Bonds used in the CME invoice-price example.
bond1 = Bond(issue_dt, Date(15, 8, 2024), 0.02375, freq, basis)
bond2 = Bond(issue_dt, Date(15, 8, 2024), 0.01875, freq, basis)

print("Futures Price       %12.6f %12.6f" % (fut_price, fut_price))

cf1 = bfut.conversion_factor(bond1)
cf2 = bfut.conversion_factor(bond2)

print("x CF                %12.4f %12.4f" % (cf1, cf2))
print("x 1000              %12.2f %12.2f" % (1000, 1000))

pip1 = bfut.principal_invoice(bond1, fut_price)
pip2 = bfut.principal_invoice(bond2, fut_price)

print("Principal invoice   %12.2f %12.2f " % (pip1, pip2))

# Convert Treasury price quotations into cash values.
price1 = 101 + 7 / 32 + 1 / 64
price2 = 98 + 1 / 32 + 0 / 64

cash1 = price1 * 1000
cash2 = price2 * 1000

print("Cash Price          %12.2f %12.2f" % (-cash1, -cash2))

# Total invoice amount includes accrued interest on the delivered bond.
tia1 = bfut.total_invoice_amount(settle_dt, bond1, fut_price)
tia2 = bfut.total_invoice_amount(settle_dt, bond2, fut_price)

print("Total Invoice price %12.2f %12.2f" % (tia1, tia2))


# ============================================================================
# 3. CME DELIVERABLE-BASKET ANALYSIS
# ============================================================================
#
# Analyse the bonds eligible for delivery into the TYZ7 Treasury futures
# contract using prices from the CME example.
#
# For each deliverable bond the example calculates:
#
#   - yield to maturity
#   - conversion factor
#   - principal invoice price
#   - total invoice amount
#   - gross basis
#   - net basis
#   - implied repo rate
#
# These measures allow the relative economics of delivering each bond
# against the futures contract to be compared.
# ============================================================================

print("\n" + "=" * 78)
print("1. BOND FUTURES CME TABLE")
print("=" * 78)

freq = FrequencyTypes.SEMI_ANNUAL
basis = DayCountTypes.ACT_ACT_ICMA
issue_dt = Date(15, 2, 2004)
yield_convention = YTMCalcType.US_TREASURY

print("TABLE 3 EXAMPLE FROM CME")
print("========================")

settle_dt = Date(10, 10, 2017)

bonds = []
prices = []
clean_prices = []

# Construct the CME deliverable basket and corresponding observed prices.

# Bond 1
bond = Bond(issue_dt, Date(15, 8, 2027), 0.0225, freq, basis)
bonds.append(bond)
prices.append(99 + 1 / 32)
clean_prices.append(99.0391)

# Bond 2
bond = Bond(issue_dt, Date(15, 5, 2027), 0.02375, freq, basis)
bonds.append(bond)
prices.append(100 + 5 / 32 + 1 / 64)
clean_prices.append(100.168)

# Bond 3
bond = Bond(issue_dt, Date(15, 2, 2027), 0.0225, freq, basis)
bonds.append(bond)
prices.append(99 + 5 / 32 + 1 / 64)
clean_prices.append(99.1641)

# Bond 4
bond = Bond(issue_dt, Date(15, 11, 2026), 0.02, freq, basis)
bonds.append(bond)
prices.append(97 + 7 / 32 + 1 / 64)
clean_prices.append(97.2305)

# Bond 5
bond = Bond(issue_dt, Date(15, 8, 2026), 0.015, freq, basis)
bonds.append(bond)
prices.append(93 + 14 / 32)
clean_prices.append(93.4414)

# Bond 6
bond = Bond(issue_dt, Date(15, 5, 2026), 0.01625, freq, basis)
bonds.append(bond)
prices.append(94 + 21 / 32 + 1 / 64)
clean_prices.append(94.6641)

# Bond 7
bond = Bond(issue_dt, Date(15, 2, 2026), 0.01625, freq, basis)
bonds.append(bond)
prices.append(94 + 29 / 32)
clean_prices.append(94.9063)

# Bond 8
bond = Bond(issue_dt, Date(15, 11, 2025), 0.0225, freq, basis)
bonds.append(bond)
prices.append(99 + 25 / 32)
clean_prices.append(99.7813)

# Bond 9
bond = Bond(issue_dt, Date(15, 8, 2025), 0.02, freq, basis)
bonds.append(bond)
prices.append(98 + 3 / 32)
clean_prices.append(98.0938)

# Bond 10
bond = Bond(issue_dt, Date(15, 5, 2025), 0.02125, freq, basis)
bonds.append(bond)
prices.append(99 + 5 / 32 + 1 / 64)
clean_prices.append(99.1719)

# Bond 11
bond = Bond(issue_dt, Date(15, 2, 2025), 0.02, freq, basis)
bonds.append(bond)
prices.append(98 + 14 / 32 + 1 / 64)
clean_prices.append(98.4531)

# Bond 12
bond = Bond(issue_dt, Date(15, 11, 2024), 0.0225, freq, basis)
bonds.append(bond)
prices.append(100 + 9 / 32 + 1 / 64)
clean_prices.append(100.3008)

# Bond 13
bond = Bond(issue_dt, Date(30, 9, 2024), 0.02125, freq, basis)
bonds.append(bond)
prices.append(99.6016)
clean_prices.append(99.6016)

# Bond 14
bond = Bond(issue_dt, Date(31, 8, 2024), 0.01875, freq, basis)
bonds.append(bond)
prices.append(98 + 1 / 32)
clean_prices.append(98.0508)

# Bond 15
bond = Bond(issue_dt, Date(15, 8, 2024), 0.02375, freq, basis)
bonds.append(bond)
prices.append(101 + 7 / 32 + 1 / 64)
clean_prices.append(101.2266)

# Bond 16
bond = Bond(issue_dt, Date(31, 7, 2024), 0.02125, freq, basis)
bonds.append(bond)
prices.append(99.6758)
clean_prices.append(99.6758)

# Bond 17
bond = Bond(issue_dt, Date(30, 6, 2024), 0.02, freq, basis)
bonds.append(bond)
prices.append(98.9336)
clean_prices.append(98.9336)

bonds.reverse()
prices.reverse()
clean_prices.reverse()

print("BOND MATURITY", "COUPON", "PRICE")
for bond, clean_price in zip(bonds, clean_prices):
    print(str(bond.maturity_dt), str(bond.cpn), clean_price)

# Recover the yield corresponding to each observed clean price.
print("BOND MATURITY", "COUPON", "YIELD")
for bond, clean_price in zip(bonds, clean_prices):
    yld = bond.yield_to_maturity(settle_dt, clean_price)
    print(str(bond.maturity_dt), str(bond.cpn), yld)

first_delivery_dt = Date(1, 12, 2017)
last_delivery_dt = Date(29, 12, 2017)

contract_size = 100000
contract_cpn = 0.06

bfut = BondFuture(
    "TYZ7",
    first_delivery_dt,
    last_delivery_dt,
    contract_size,
    contract_cpn,
)

# Conversion factors standardise the bonds in the delivery basket against
# the notional coupon of the futures contract.
print("BOND MATURITY", "COUPON", "CF")
for bond in bonds:
    cf = bfut.conversion_factor(bond)
    print(str(bond.maturity_dt), str(bond.cpn), cf)

futures_price = 125.265625

# Principal invoice price before accrued interest.
print("BOND MATURITY", "PRINCIPAL INVOICE PRICE")
for bond in bonds:
    pip = bfut.principal_invoice(bond, futures_price)
    print(str(bond.maturity_dt), pip)

# Total invoice amount includes the accrued interest payable on delivery.
print("BOND MATURITY", "TOTAL INVOICE AMOUNT")
for bond in bonds:
    tia = bfut.total_invoice_amount(settle_dt, bond, futures_price)
    print(str(bond.maturity_dt), tia)

# Implied repo measures the return from the cash-and-carry trade in which
# the bond is purchased, financed and subsequently delivered into futures.
print("BOND MATURITY", "IMPLIED REPO RATE")
for bond, clean_price in zip(bonds, clean_prices):
    repo_rate = bfut.implied_repo_rate(
        bond,
        settle_dt,
        clean_price,
        futures_price,
    )
    print(str(bond.maturity_dt), repo_rate)

# Identify the cheapest-to-deliver bond from the delivery basket.
ctd = bfut.ctd(bonds, prices, futures_price)

print("CTD MATURITY", "CTD COUPON")
print(str(ctd.maturity_dt), ctd.cpn)

results = []

# Financing rate used when calculating the net basis.
repo_rate = 0.015

for bond, clean_price in zip(bonds, clean_prices):

    cf = bfut.conversion_factor(bond)

    yld = bond.yield_to_maturity(
        settle_dt,
        clean_price,
        yield_convention,
    )

    del_years = bfut.delivery_years(bond)

    pip = bfut.principal_invoice(
        bond,
        futures_price,
    )

    tia = bfut.total_invoice_amount(
        settle_dt,
        bond,
        futures_price,
    )

    gross_basis = bfut.gross_basis(
        bond,
        clean_price,
        futures_price,
    )

    net_basis = bfut.net_basis(
        bond,
        settle_dt,
        clean_price,
        futures_price,
        repo_rate,
    )

    irr = bfut.implied_repo_rate(
        bond,
        settle_dt,
        clean_price,
        futures_price,
    )

    # Placeholder retained from the original example.
    fut_dv01 = 0.0

    results.append(
        {
            "Coupon": bond.cpn * 100,
            "Maturity": bond.maturity_dt,
            "Clean": clean_price,
            "DelYrs": del_years,
            "TCF": cf,
            "PIP": pip,
            "TIA": tia,
            "GROSS_BASIS": gross_basis,
            "NET_BASIS": net_basis,
            "FUT_DV01": fut_dv01,
            "Yield (%)": yld * 100,
            "IRR (%)": irr * 100,
        }
    )

# Present the delivery basket ranked by implied repo rate.
df = pd.DataFrame(results)
df = df.sort_values(
    by="IRR (%)",
    ascending=False,
).reset_index(drop=True)

formatted_df = df.copy()

formatted_df["Coupon"] = formatted_df["Coupon"].map("{:.3f}".format)
formatted_df["DelYrs"] = formatted_df["DelYrs"].map("{:.2f}".format)
formatted_df["TCF"] = formatted_df["TCF"].map("{:.4f}".format)
formatted_df["PIP"] = formatted_df["PIP"].map("{:.2f}".format)
formatted_df["TIA"] = formatted_df["TIA"].map("{:.2f}".format)
formatted_df["GROSS_BASIS"] = formatted_df["GROSS_BASIS"].map("{:.3f}".format)
formatted_df["NET_BASIS"] = formatted_df["NET_BASIS"].map("{:.3f}".format)
formatted_df["FUT_DV01"] = formatted_df["FUT_DV01"].map("{:.3f}".format)
formatted_df["IRR (%)"] = formatted_df["IRR (%)"].map("{:.3f}".format)
formatted_df["Yield (%)"] = formatted_df["Yield (%)"].map("{:.3f}".format)

print(formatted_df.to_string(index=False))


# ============================================================================
# 4. BLOOMBERG DELIVERABLE-BASKET ANALYSIS
# ============================================================================
#
# Repeat the TYZ7 deliverable-basket analysis using the Bloomberg prices
# supplied in the original example.
#
# The same bond-futures measures are calculated so that the observed cash
# prices can be compared on a consistent basis:
#
#   - yield to maturity
#   - conversion factor
#   - principal invoice price
#   - total invoice amount
#   - gross basis
#   - net basis
#   - implied repo rate
#   - cheapest-to-deliver bond
# ============================================================================

print("\n" + "=" * 78)
print("2. BOND FUTURES BBG TABLE")
print("=" * 78)

freq = FrequencyTypes.SEMI_ANNUAL
basis = DayCountTypes.ACT_ACT_ICMA
issue_dt = Date(15, 2, 2004)
yield_convention = YTMCalcType.US_TREASURY

print("TABLE 3 EXAMPLE FROM CME")
print("========================")

settle_dt = Date(10, 10, 2017)

bonds = []
prices = []
new_prices = []

# Bloomberg clean prices for the deliverable bonds.

bond = Bond(issue_dt, Date(31, 7, 2024), 0.02125, freq, basis)
bonds.append(bond)
new_prices.append(98.5703)

bond = Bond(issue_dt, Date(30, 6, 2024), 0.020, freq, basis)
bonds.append(bond)
new_prices.append(97.8516)

bond = Bond(issue_dt, Date(31, 8, 2024), 0.01875, freq, basis)
bonds.append(bond)
new_prices.append(97.0469)

bond = Bond(issue_dt, Date(30, 11, 2024), 0.02125, freq, basis)
bonds.append(bond)
new_prices.append(98.4297)

bond = Bond(issue_dt, Date(31, 10, 2024), 0.0225, freq, basis)
bonds.append(bond)
new_prices.append(99.2578)

bond = Bond(issue_dt, Date(15, 11, 2024), 0.0225, freq, basis)
bonds.append(bond)
new_prices.append(99.2188)

bond = Bond(issue_dt, Date(30, 9, 2024), 0.02125, freq, basis)
bonds.append(bond)
new_prices.append(98.4844)

bond = Bond(issue_dt, Date(15, 8, 2024), 0.02375, freq, basis)
bonds.append(bond)
new_prices.append(101.2266)

print("BOND MATURITY", "COUPON", "PRICE")
for bond, clean_price in zip(bonds, new_prices):
    print(str(bond.maturity_dt), str(bond.cpn), clean_price)

# Convert each observed clean price into its corresponding Treasury yield.
print("BOND MATURITY", "COUPON", "YIELD")
for bond, clean_price in zip(bonds, new_prices):
    yld = bond.yield_to_maturity(
        settle_dt,
        clean_price,
        yield_convention,
    )
    print(str(bond.maturity_dt), str(bond.cpn), yld)

first_delivery_dt = Date(1, 12, 2017)
last_delivery_dt = Date(29, 12, 2017)

contract_size = 100000
contract_cpn = 0.06

bfut = BondFuture(
    "TYZ7",
    first_delivery_dt,
    last_delivery_dt,
    contract_size,
    contract_cpn,
)

# Calculate the conversion factor for each deliverable bond.
print("BOND MATURITY", "COUPON", "CF")
for bond in bonds:
    cf = bfut.conversion_factor(bond)
    print(str(bond.maturity_dt), str(bond.cpn), cf)

futures_price = 125.265625

# Calculate principal invoice prices at the observed futures price.
print("BOND MATURITY", "PRINCIPAL INVOICE PRICE")
for bond in bonds:
    pip = bfut.principal_invoice(
        bond,
        futures_price,
    )
    print(str(bond.maturity_dt), pip)

# Add accrued interest to obtain the total invoice amount.
print("BOND MATURITY", "TOTAL INVOICE AMOUNT")
for bond in bonds:
    tia = bfut.total_invoice_amount(
        settle_dt,
        bond,
        futures_price,
    )
    print(str(bond.maturity_dt), tia)

# Determine the cheapest-to-deliver bond using the observed cash prices.
ctd = bfut.ctd(
    bonds,
    new_prices,
    futures_price,
)

print("CTD MATURITY", "CTD COUPON")
print(str(ctd.maturity_dt), ctd.cpn)

results = []

# Repo rate used in the net-basis calculation.
repo_rate = 0.01499

for bond, clean_price in zip(bonds, new_prices):

    mid_price = clean_price

    cf = bfut.conversion_factor(bond)

    yld = bond.yield_to_maturity(
        settle_dt,
        mid_price,
        yield_convention,
    )

    del_years = bfut.delivery_years(bond)

    pip = bfut.principal_invoice(
        bond,
        futures_price,
    )

    tia = bfut.total_invoice_amount(
        settle_dt,
        bond,
        futures_price,
    )

    gross_basis = bfut.gross_basis(
        bond,
        mid_price,
        futures_price,
    )

    irr = bfut.implied_repo_rate(
        bond,
        settle_dt,
        mid_price,
        futures_price,
    )

    net_basis = bfut.net_basis(
        bond,
        settle_dt,
        clean_price,
        futures_price,
        repo_rate,
    )

    results.append(
        {
            "Coupon": bond.cpn * 100,
            "Maturity": bond.maturity_dt,
            "Mid": mid_price,
            "DelYrs": del_years,
            "TCF": cf,
            "PIP": pip,
            "TIA": tia,
            "GROSS_BASIS": gross_basis,
            "NET_BASIS": net_basis,
            "Yield (%)": yld * 100,
            "IRR (%)": irr * 100,
        }
    )

# Format the final comparison table without modifying the underlying values.
df = pd.DataFrame(results)
formatted_df = df.copy()

formatted_df["Coupon"] = formatted_df["Coupon"].map("{:.3f}".format)
formatted_df["DelYrs"] = formatted_df["DelYrs"].map("{:.2f}".format)
formatted_df["TCF"] = formatted_df["TCF"].map("{:.4f}".format)
formatted_df["PIP"] = formatted_df["PIP"].map("{:.2f}".format)
formatted_df["TIA"] = formatted_df["TIA"].map("{:.2f}".format)
formatted_df["GROSS_BASIS"] = formatted_df["GROSS_BASIS"].map("{:.3f}".format)
formatted_df["NET_BASIS"] = formatted_df["NET_BASIS"].map("{:.3f}".format)
formatted_df["IRR (%)"] = formatted_df["IRR (%)"].map("{:.3f}".format)
formatted_df["Yield (%)"] = formatted_df["Yield (%)"].map("{:.6f}".format)

print(formatted_df.to_string(index=False))
