# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import time
import numpy as np

import matplotlib.pyplot as plt

from financepy.models.heston import Heston, HestonValueTypes
from financepy.utils.global_types import OptionTypes, HestonNumericalSchemeTypes
from financepy.utils.global_vars import G_DAYS_IN_YEAR
from financepy.utils.date import Date

# ============================================================================
# FINANCEPY EXAMPLES - Heston
# ============================================================================



PLOT = False

########################################################################################




########################################################################################




########################################################################################




########################################################################################

# ============================================================================
# 1. ANALYTICAL MODELS
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. ANALYTICAL MODELS")
print("=" * 78)

value_dt = Date(1, 1, 2015)
expiry_dt = Date(1, 4, 2015)
v0 = 0.05  # initial variance of volatility
theta = 0.05  # long term variance
kappa = 2.0  # speed of variance reversion
sigma = 0.10  # volatility of variance
rho = -0.9  # correlation
interest_rate = 0.05
dividend_yield = 0.01
seed = 2838

tau = (expiry_dt - value_dt) / G_DAYS_IN_YEAR
opt_type = OptionTypes.EUROPEAN_CALL.value

num_steps = 100
num_paths = 20000
stock_price = 100.0

print(
    "TIME",
    "RHO",
    "SIGMA",
    "K",
    "MC",
    "GATH",
    "LEWROU",
    "LEWIS",
    "WEBER",
    "MCERR",
)

for sigma in [0.5, 0.75, 1.0]:
    for rho in [-0.9, -0.5, 0.0]:
        heston_model = Heston(v0, kappa, theta, sigma, rho)

        for strike_price in np.linspace(95, 105, 3):
            value_mc_heston = heston_model.value_mc(
                stock_price,
                tau,
                strike_price,
                opt_type,
                interest_rate,
                dividend_yield,
                num_paths,
                num_steps,
                seed,
            )

            start = time.time()
            value_gatheral = heston_model.value(
                stock_price,
                tau,
                strike_price,
                opt_type,
                interest_rate,
                dividend_yield,
                method=HestonValueTypes.GATHERAL,
            )

            value_lewis_rouah = heston_model.value(
                stock_price,
                tau,
                strike_price,
                opt_type,
                interest_rate,
                dividend_yield,
                method=HestonValueTypes.LEWIS_ROUAH,
            )

            value_lewis = heston_model.value(
                stock_price,
                tau,
                strike_price,
                opt_type,
                interest_rate,
                dividend_yield,
                method=HestonValueTypes.LEWIS,
            )

            value_weber = heston_model.value(
                stock_price,
                tau,
                strike_price,
                opt_type,
                interest_rate,
                dividend_yield,
                method=HestonValueTypes.WEBER,
            )

    err = value_mc_heston - value_weber
    end = time.time()
    elapsed = end - start
    print(
        f"{elapsed:6.3f}",
        f"{rho: 7.5f}",
        f"{sigma:7.5f}",
        f"{strike_price:7.2f}",
        f"{value_mc_heston:12.9f}",
        f"{value_gatheral:12.9f}",
        f"{value_lewis_rouah:12.9f}",
        f"{value_lewis:12.9f}",
        f"{value_weber:12.9f}",
        f"{err:12.9f}",
    )

# ============================================================================
# 2. MONTE CARLO
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("2. MONTE CARLO")
print("=" * 78)

value_dt = Date(1, 1, 2015)
expiry_dt = Date(1, 1, 2016)
v0 = 0.04  # initial variance of volatility
theta = 0.04  # long term variance
kappa = 2.0  # speed of variance reversion
sigma = 1.0  # volatility of variance
rho = -0.9  # correlation
interest_rate = 0.05
dividend_yield = 0.01
seed = 238

stock_price = 100.0

print(
    "TIME",
    "RHO",
    "SIGMA",
    "K",
    "NSTEPS",
    "NPATHS",
    "FORMULA",
    "EULER_ERR",
    "EULLOG_ERR",
    "QE_ERR",
)

tau = (expiry_dt - value_dt) / G_DAYS_IN_YEAR
opt_type = OptionTypes.EUROPEAN_CALL.value

for strike_price in np.linspace(95, 105, 3):
    for num_steps in [25, 50]:
        for num_paths in [10000, 20000]:

            heston_model = Heston(v0, kappa, theta, sigma, rho)

            value_weber = heston_model.value(
                stock_price,
                tau,
                strike_price,
                opt_type,
                interest_rate,
                dividend_yield,
                method=HestonValueTypes.WEBER,
            )

            start = time.time()

            value_mc_euler = heston_model.value_mc(
                stock_price,
                tau,
                strike_price,
                opt_type,
                interest_rate,
                dividend_yield,
                num_paths,
                num_steps,
                seed,
                HestonNumericalSchemeTypes.EULER,
            )
            value_mc_euler_log = heston_model.value_mc(
                stock_price,
                tau,
                strike_price,
                opt_type,
                interest_rate,
                dividend_yield,
                num_paths,
                num_steps,
                seed,
                HestonNumericalSchemeTypes.EULERLOG,
            )
            value_mc_quadexp = heston_model.value_mc(
                stock_price,
                tau,
                strike_price,
                opt_type,
                interest_rate,
                dividend_yield,
                num_paths,
                num_steps,
                seed,
                HestonNumericalSchemeTypes.QUADEXP,
            )

            err_euler = value_mc_euler - value_weber
            err_euler_log = value_mc_euler_log - value_weber
            err_quadexp = value_mc_quadexp - value_weber

            end = time.time()
            elapsed = end - start

            print(
                elapsed,
                rho,
                sigma,
                strike_price,
                num_steps,
                num_paths,
                value_weber,
                err_euler,
                err_euler_log,
                err_quadexp,
            )

# ============================================================================
# 3. HESTON VOLATILITY SMILE
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("3. HESTON VOLATILITY SMILE")
print("=" * 78)

model = Heston(
    v0=0.04,
    kappa=2.0,
    theta=0.04,
    xi=0.50,
    rho=-0.70,
)

model = Heston(
    v0=0.04,
    kappa=2.0,
    theta=0.04,
    xi=0.80,
    rho=0.0,
)

t_exp = 1.0
stock_price = 100.0
interest_rate = 0.05
dividend_yield = 0.02

strikes = np.linspace(60.0, 140.0, 41)

vols = model.volatility_smile(
    t_exp,
    strikes,
    stock_price,
    interest_rate,
    dividend_yield,
    HestonValueTypes.LEWIS,
)

# Only check finite values after diagnostics have been printed
valid = np.isfinite(vols)

assert np.all(valid)
assert np.all(vols[valid] > 0.0)

# Plot
if PLOT == True:
    plt.figure()
    plt.plot(
        strikes[valid],
        100.0 * vols[valid],
        marker="o",
        markersize=3,
    )

    plt.axvline(stock_price, linestyle="--")
    plt.xlabel("Strike")
    plt.ylabel("Implied Volatility (%)")
    plt.title("Heston Implied Volatility Smile")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

