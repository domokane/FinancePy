# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import time
import numpy as np


from financepy.models.cir_montecarlo import zero_price_mc, zero_price
from financepy.utils.global_types import CIRNumericalSchemeTypes

# ============================================================================
# FINANCEPY EXAMPLES - Model Cir
# ============================================================================


########################################################################################




########################################################################################

# ============================================================================
# 1. FIN MODEL RATES CIR
# ============================================================================
# What this section demonstrates:
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. FIN MODEL RATES CIR")
print("=" * 78)

r0 = 0.05
a = 0.20
b = 0.05
sigma = 0.20
t = 5.0

num_paths = 2000
dt = 0.05
seed = 1968

print(
    "MATURITY",
    "TIME",
    "FORMULA",
    "EULER",
    "LOGNORM",
    "MILSTEIN",
    "KJ",
    "EXACT",
)

for t in np.linspace(0, 10, 21):

    start = time.time()
    p = zero_price(r0, a, b, sigma, t)
    p_mc1 = zero_price_mc(
        r0,
        a,
        b,
        sigma,
        t,
        dt,
        num_paths,
        seed,
        CIRNumericalSchemeTypes.EULER.value,
    )
    p_mc2 = zero_price_mc(
        r0,
        a,
        b,
        sigma,
        t,
        dt,
        num_paths,
        seed,
        CIRNumericalSchemeTypes.LOGNORMAL.value,
    )
    p_mc3 = zero_price_mc(
        r0,
        a,
        b,
        sigma,
        t,
        dt,
        num_paths,
        seed,
        CIRNumericalSchemeTypes.MILSTEIN.value,
    )
    p_mc4 = zero_price_mc(
        r0,
        a,
        b,
        sigma,
        t,
        dt,
        num_paths,
        seed,
        CIRNumericalSchemeTypes.KAHLJACKEL.value,
    )
    p_mc5 = zero_price_mc(
        r0,
        a,
        b,
        sigma,
        t,
        dt,
        num_paths,
        seed,
        CIRNumericalSchemeTypes.EXACT.value,
    )
    end = time.time()
    elapsed = end - start
    print(t, elapsed, p, p_mc1, p_mc2, p_mc3, p_mc4, p_mc5)

