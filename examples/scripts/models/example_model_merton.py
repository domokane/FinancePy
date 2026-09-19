# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import numpy as np

from financepy.models.merton_firm_mkt import MertonFirmMkt
from financepy.models.merton_firm import MertonFirm
import matplotlib.pyplot as plt

# ============================================================================
# FINANCEPY EXAMPLES - MertonFirmMkt
# ============================================================================


########################################################################################




########################################################################################

# ============================================================================
# 1. FIN MODEL MERTON CREDIT
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("1. FIN MODEL MERTON CREDIT")
print("=" * 78)

equity_value = np.array([2.6406, 2.6817, 3.9770, 2.9470, 2.5280])
equity_vol = np.array([0.7103, 0.3929, 0.3121, 0.4595, 0.6181])
bond_face = np.array([4.0, 3.5, 3.5, 3.2, 4.0])
risk_free_rate = 0.05
asset_growth_rate = np.array([0.0306, 0.0300, 0.0310, 0.0302, 0.0305])
years_to_maturity = 1.0

model_mkt = MertonFirmMkt(
    equity_value,
    bond_face,
    years_to_maturity,
    risk_free_rate,
    asset_growth_rate,
    equity_vol,
)

print("MERTON MARKET MODEL", "VALUE")

print("ASSET VALUE", model_mkt.asset_value())
print("EQUITY VALUE", model_mkt.equity_value())
print("DEBT VALUE", model_mkt.debt_value())

print("ASSET VOLATILITY", model_mkt.asset_volatility())
print("EQUITY VOL", model_mkt.equity_volatility())

print("CREDIT SPREAD", model_mkt.credit_spread())
print("ASSET TO DEBT", model_mkt.asset_to_debt_ratio())
print("RISK NEUTRAL PROB DEFAULT", model_mkt.risk_neutral_default_probability())
print("PHYSICAL PROB DEFAULT", model_mkt.physical_default_probability())
print("DISTANCE DEFAULT", model_mkt.distance_to_default())

# -------------------------------------------------------------------------
# Check that the inferred asset quantities reproduce the market inputs.
# -------------------------------------------------------------------------

assert np.allclose(
    model_mkt.equity_value(),
    equity_value,
    rtol=1.0e-8,
    atol=1.0e-10,
)

assert np.allclose(
    model_mkt.equity_volatility(),
    equity_vol,
    rtol=1.0e-8,
    atol=1.0e-10,
)

# -------------------------------------------------------------------------
# Pass inferred A and sigma_A into the basic Merton model.
# -------------------------------------------------------------------------

asset_value = model_mkt.asset_value()
asset_vol = model_mkt.asset_volatility()

model = MertonFirm(
    asset_value,
    bond_face,
    years_to_maturity,
    risk_free_rate,
    asset_growth_rate,
    asset_vol,
)

print("BASIC MERTON MODEL", "VALUE")

print("ASSET VALUE", model.asset_value())
print("EQUITY VALUE", model.equity_value())
print("DEBT VALUE", model.debt_value())

print("ASSET VOLATILITY", model.asset_volatility())
print("EQUITY VOL", model.equity_volatility())

print("CREDIT SPREAD", model.credit_spread())
print("ASSET TO DEBT", model.asset_to_debt_ratio())
print("RISK NEUTRAL DEFAULT PROB", model.risk_neutral_default_probability())
print("PHYSICAL DEFAULT PROB", model.physical_default_probability())
print("DISTANCE DEFAULT", model.distance_to_default())

# -------------------------------------------------------------------------
# Basic Merton balance-sheet identity:
#
# A = E + D
# -------------------------------------------------------------------------

assert np.allclose(
    model.asset_value(),
    model.equity_value() + model.debt_value(),
    rtol=1.0e-12,
    atol=1.0e-12,
)

# Market and basic implementations should agree after inversion.

assert np.allclose(
    model.equity_value(),
    model_mkt.equity_value(),
    rtol=1.0e-10,
    atol=1.0e-12,
)

assert np.allclose(
    model.debt_value(),
    model_mkt.debt_value(),
    rtol=1.0e-10,
    atol=1.0e-12,
)

# -------------------------------------------------------------------------
# Scalar example.
# -------------------------------------------------------------------------

asset_value = 140.0
bond_face = 100.0
years_to_maturity = 1.0
risk_free_rate = 0.05
asset_growth_rate = 0.05
asset_vol = 0.20

model = MertonFirm(
    asset_value,
    bond_face,
    years_to_maturity,
    risk_free_rate,
    asset_growth_rate,
    asset_vol,
)

print("BASIC MERTON MODEL SCALAR", "VALUE")

print("ASSET VALUE", model.asset_value())
print("EQUITY VALUE", model.equity_value())
print("DEBT VALUE", model.debt_value())

print("ASSET VOLATILITY", model.asset_volatility())
print("EQUITY VOLATILITY", model.equity_volatility())

print("CREDIT SPREAD", model.credit_spread())
print("ASSET TO DEBT", model.asset_to_debt_ratio())
print("RISK NEUTRAL DEFAULT PROB", model.risk_neutral_default_probability())
print("PHYSICAL DEFAULT PROB", model.physical_default_probability())
print("DISTANCE DEFAULT", model.distance_to_default())

assert np.allclose(
    model.asset_value(),
    model.equity_value() + model.debt_value(),
    rtol=1.0e-12,
    atol=1.0e-12,
)

# =============================================================================
# 2. VISUALISE DEFAULT RISK ACROSS THE MARKET EXAMPLES
# =============================================================================
# Distance to default and default probability summarize the same structural
# credit story in different units. The chart makes cross-company differences
# visible without relying on a dense vector printout.
plot_labels = [f"Firm {i + 1}" for i in range(len(equity_value))]
plot_pd = np.asarray(model_mkt.risk_neutral_default_probability()) * 100.0
plt.figure()
plt.bar(plot_labels, plot_pd)
plt.xlabel("Firm")
plt.ylabel("Risk-neutral default probability (%)")
plt.title("Merton model: risk-neutral default probability by firm")
plt.grid(True, axis="y")
plt.tight_layout()
plt.show()
