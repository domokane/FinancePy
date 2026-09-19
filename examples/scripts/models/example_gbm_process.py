# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import numpy as np


from financepy.models.gbm_process_simulator import get_assets_paths_times
from financepy.utils.math import corr_matrix_generator
import matplotlib.pyplot as plt

# ============================================================================
# FINANCEPY EXAMPLES - Gbm Process
# ============================================================================


########################################################################################




########################################################################################

# ============================================================================
# 1. FIN GBM PROCESS
# ============================================================================
# What this section demonstrates:
# Runs the original FinancePy calculation with explicit inputs so the numerical result and the effect of the chosen assumptions can be inspected.

print("\n" + "=" * 78)
print("1. FIN GBM PROCESS")
print("=" * 78)

num_assets = 3
num_paths = 6
num_time_steps = 1
t = 1.0
mus = 0.03 * np.ones(num_assets)
stock_prices = 100.0 * np.ones(num_assets)
volatilities = 0.2 * np.ones(num_assets)
rho = 0.8
corr_matrix = corr_matrix_generator(rho, num_assets)
seed = 1912

times, paths = get_assets_paths_times(
    num_assets,
    num_paths,
    num_time_steps,
    t,
    mus,
    stock_prices,
    volatilities,
    corr_matrix,
    seed,
)

# =============================================================================
# 2. VISUALISE SIMULATED CORRELATED ASSET OUTCOMES
# =============================================================================
# With one time step, the most useful view is the distribution of terminal
# values across assets and paths. Assets share a high positive correlation, so
# their simulated outcomes tend to move in the same broad direction.
plot_paths = np.asarray(paths)
plt.figure()
for asset_index in range(num_assets):
    terminal_values = plot_paths[:, asset_index, -1] if plot_paths.ndim == 3 else plot_paths[asset_index]
    plt.plot(range(1, len(terminal_values) + 1), terminal_values,
             marker="o", label=f"Asset {asset_index + 1}")
plt.xlabel("Simulation path")
plt.ylabel("Terminal asset value")
plt.title("Correlated GBM terminal values")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
