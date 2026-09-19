# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

# Allow this example to run directly from its category folder.

from financepy.utils.currency import CurrencyTypes
from financepy.utils.amount import Amount

# ============================================================================
# FINANCEPY EXAMPLES - Amount
# ============================================================================



########################################################################################




########################################################################################

# ============================================================================
# 1. AMOUNT
# ============================================================================
# What this section demonstrates:
# Runs the original FinancePy calculation with explicit inputs so the numerical result and the effect of the chosen assumptions can be inspected.

print("\n" + "=" * 78)
print("1. AMOUNT")
print("=" * 78)

print("LABEL", "AMOUNT")
x = Amount(101000.232, CurrencyTypes.USD)

print("Amount", x)

x = Amount(101000.232, CurrencyTypes.CAD)

print("Amount", x)

