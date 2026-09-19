# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import numpy as np
import matplotlib.pyplot as plt


from financepy.models.sabr import SABR
from financepy.models.sabr_shifted import SABRShifted

# ============================================================================
# FINANCEPY EXAMPLES - Model Black Sabr Hw
# ============================================================================



PLOT_GRAPHS = False

########################################################################################




########################################################################################




########################################################################################




########################################################################################

# ============================================================================
# 1. FIN SABR
# ============================================================================
# What this section demonstrates:
# Runs the original FinancePy calculation with explicit inputs so the numerical result and the effect of the chosen assumptions can be inspected.

print("\n" + "=" * 78)
print("1. FIN SABR")
print("=" * 78)

strikes = np.linspace(0.01, 0.06, 10)

alpha = 0.060277
beta = 0.5
rho = 0.2097
nu = 0.75091
model1 = SABR(alpha, beta, rho, nu)

alpha = 0.058484
beta = 0.5
rho = 0.20568
nu = 0.79647
model2 = SABR(alpha, beta, rho, nu)

f = 0.0350
t = 1.0

vols1 = model1.black_vol(f, strikes, t)
vols2 = model2.black_vol(f, strikes, t)

if PLOT_GRAPHS:
    plt.figure()
    plt.plot(strikes, vols1)
    plt.plot(strikes, vols2)
    plt.title("SABR")

# ============================================================================
# 2. FIN SHIFTED SABR SIMPLE
# ============================================================================
# What this section demonstrates:
# Runs the original FinancePy calculation with explicit inputs so the numerical result and the effect of the chosen assumptions can be inspected.

print("\n" + "=" * 78)
print("2. FIN SHIFTED SABR SIMPLE")
print("=" * 78)

strikes = np.linspace(0.01, 0.06, 10)

alpha = 0.060277
beta = 0.5
rho = 0.2097
nu = 0.75091
model1 = SABRShifted(alpha, beta, rho, nu, 0.0)

alpha = 0.058484
beta = 0.5
rho = 0.20568
nu = 0.79647
model2 = SABRShifted(alpha, beta, rho, nu, 0.0)

f = 0.0350
t = 1.0

vols1 = model1.black_vol(f, strikes, t)
vols2 = model2.black_vol(f, strikes, t)

if PLOT_GRAPHS:
    plt.figure()
    plt.plot(strikes, vols1)
    plt.plot(strikes, vols2)
    plt.title("Shifted SIMPLE SABR")

# ============================================================================
# 3. FIN SHIFTED SABR
# ============================================================================
# What this section demonstrates:
# Runs the original FinancePy calculation with explicit inputs so the numerical result and the effect of the chosen assumptions can be inspected.

print("\n" + "=" * 78)
print("3. FIN SHIFTED SABR")
print("=" * 78)

strikes = np.linspace(-0.006, 0.016, 10)

alpha = 0.013345
beta = 0.5
rho = 0.46698
nu = 0.49861
shift = 0.008

model = SABRShifted(alpha, beta, rho, nu, shift)

f = 0.0006384
t = 1.0

vols = model.black_vol(f, strikes, t)

if PLOT_GRAPHS:
    plt.figure()
    plt.plot(strikes, vols)
    plt.title("SHIFTED SABR")

