# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.


import time
import numpy as np


from financepy.utils.date import Date
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.products.bonds.bond import Bond
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.global_vars import G_DAYS_IN_YEAR
from financepy.utils.helpers import print_tree
from financepy.models.bk_tree import BKTree
from financepy.utils.global_types import ExerciseTypes

# ============================================================================
# FINANCEPY EXAMPLES - Model Rates Bk
# ============================================================================

########################################################################################




########################################################################################




########################################################################################

# ============================================================================
# 1. BK EXAMPLE ONE
# ============================================================================
# What this section demonstrates:
# Runs the original FinancePy calculation with explicit inputs so the numerical result and the effect of the chosen assumptions can be inspected.

print("\n" + "=" * 78)
print("1. BK EXAMPLE ONE")
print("=" * 78)

print("=== HULL INITIAL EXAMPLE SECTION 28.7 ED 6 PG 668 ====")

times = [0.0, 0.5000, 1.00000, 1.50000, 2.00000, 2.500000, 3.00000]
zeros = [0.03, 0.0343, 0.03824, 0.04183, 0.04512, 0.048512, 0.05086]
times = np.array(times)
zeros = np.array(zeros)
dfs = np.exp(-zeros * times)

start_dt = Date(1, 12, 2019)
end_dt = Date(1, 6, 2021)
sigma = 0.25
a = 0.22
num_time_steps = 3
t_mat = (end_dt - start_dt) / G_DAYS_IN_YEAR
model = BKTree(sigma, a, num_time_steps)
model.build_tree(t_mat, times, dfs)

# Agrees with Figure 28.10 - Not exact as we have dt not exactly 0.50
if num_time_steps < 5:
    print("LABEL", "VALUE")
    print("QTREE", model.qq)
    print("RTREE", model.rt)
    #        print_tree(model._rt)
    print("PU AT LAST TIME", model.pu)
    print("PDM AT LAST TIME", model.pm)
    print("PD AT LAST TIME", model.pd)

# ============================================================================
# 2. BK EXAMPLE TWO
# ============================================================================
# What this section demonstrates:
# Values cash flows from the curve and reports the clean quoted price.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("2. BK EXAMPLE TWO")
print("=" * 78)

settle_dt = Date(1, 12, 2019)
issue_dt = Date(1, 12, 2018)
expiry_dt = settle_dt.add_tenor("18m")
maturity_dt = settle_dt.add_tenor("10Y")
coupon = 0.05
freq_type = FrequencyTypes.SEMI_ANNUAL
dc_type = DayCountTypes.ACT_ACT_ICMA
bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

cpn_times = []
cpn_flows = []
cpn = bond.cpn / bond.freq
num_flows = len(bond.cpn_dts)

for i in range(1, num_flows):
    pcd = bond.cpn_dts[i - 1]
    ncd = bond.cpn_dts[i]
    if pcd < settle_dt and ncd > settle_dt:
        flow_time = (pcd - settle_dt) / G_DAYS_IN_YEAR
        cpn_times.append(flow_time)
        cpn_flows.append(cpn)

for flow_dt in bond.cpn_dts:
    if flow_dt > settle_dt:
        flow_time = (flow_dt - settle_dt) / G_DAYS_IN_YEAR
        cpn_times.append(flow_time)
        cpn_flows.append(cpn)

cpn_times = np.array(cpn_times)
cpn_flows = np.array(cpn_flows)

strike_price = 105.0
face = 100.0

t_mat = (maturity_dt - settle_dt) / G_DAYS_IN_YEAR
t_exp = (expiry_dt - settle_dt) / G_DAYS_IN_YEAR
times = np.linspace(0, t_mat, 11)
dates = settle_dt.add_years(times)
dfs = np.exp(-0.05 * times)
curve = DiscountCurve(settle_dt, dates, dfs)

price = bond.clean_price_from_discount_curve(settle_dt, curve)
print("LABEL", "VALUE")
print("Fixed Income Price:", price)

sigma = 0.20
a = 0.05
num_time_steps = 26

model = BKTree(sigma, a, num_time_steps)
model.build_tree(t_mat, times, dfs)
exercise_type = ExerciseTypes.AMERICAN
v = model.bond_option(t_exp, strike_price, face, cpn_times, cpn_flows, exercise_type)

# Test convergence
num_steps_list = [100, 200, 300, 500, 1000]
exercise_type = ExerciseTypes.AMERICAN

print("time_steps", "TIME", "VALUE")
tree_vector = []
for num_time_steps in num_steps_list:
    start = time.time()
    model = BKTree(sigma, a, num_time_steps)
    model.build_tree(t_mat, times, dfs)
    v = model.bond_option(t_exp, strike_price, face, cpn_times, cpn_flows, exercise_type)
    end = time.time()
    period = end - start
    tree_vector.append(v)
    print(num_time_steps, period, v)

#    plt.plot(num_steps_list, tree_vector)

# Value in Hill converges to 0.699 with 100 time steps while I get 0.700

if False:
    print("RT")
    print_tree(model.rt, 5)
    print("Q")
    print_tree(model.qq, 5)

