# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import time

import matplotlib.pyplot as plt
import numpy as np


from financepy.models.bdt_tree import BDTTree
from financepy.utils.global_types import OptionTypes
from financepy.products.bonds.bond_option import BondOption
from financepy.utils.global_vars import G_DAYS_IN_YEAR
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.bonds.bond import Bond
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.utils.date import Date

# ============================================================================
# FINANCEPY EXAMPLES - BondOption
# ============================================================================




PLOT_GRAPHS = False

########################################################################################




########################################################################################


def test_bond_option_american_convergence_one():

    # Build discount curve
    settle_dt = Date(1, 12, 2019)
    discount_curve = FlatDiscountCurve(settle_dt, 0.05)

    # Bond details
    issue_dt = Date(1, 9, 2010)
    maturity_dt = Date(1, 9, 2025)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    # Option Details
    expiry_dt = Date(1, 12, 2020)
    strike_price = 100.0

    print("TIME", "N", "PUT_AMER", "PUT_EUR", "CALL_AME", "CALL_EUR")

    time_steps = range(30, 100, 1)

    for num_time_steps in time_steps:

        sigma = 0.20

        start = time.time()

        opt_type = OptionTypes.AMERICAN_PUT
        bond_option1 = BondOption(bond, expiry_dt, strike_price, opt_type)
        model = BDTTree(sigma, num_time_steps)
        v1put = bond_option1.value(settle_dt, discount_curve, model)

        opt_type = OptionTypes.EUROPEAN_PUT
        bond_option2 = BondOption(bond, expiry_dt, strike_price, opt_type)
        model = BDTTree(sigma, num_time_steps)
        v2put = bond_option2.value(settle_dt, discount_curve, model)

        opt_type = OptionTypes.AMERICAN_CALL
        bond_option1 = BondOption(bond, expiry_dt, strike_price, opt_type)
        model = BDTTree(sigma, num_time_steps)
        v1call = bond_option1.value(settle_dt, discount_curve, model)

        opt_type = OptionTypes.EUROPEAN_CALL
        bond_option2 = BondOption(bond, expiry_dt, strike_price, opt_type)
        model = BDTTree(sigma, num_time_steps)
        v2call = bond_option2.value(settle_dt, discount_curve, model)

        end = time.time()

        period = end - start

        print(period, num_time_steps, v1put, v2put, v1call, v2call)


########################################################################################




########################################################################################




########################################################################################

# test_bond_option_american_convergence_one()

# ============================================================================
# 1. BOND OPTION ZEROVOL CONVERGENCE
# ============================================================================
# What this section demonstrates:
# Values cash flows from the curve and reports the clean quoted price.
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. BOND OPTION ZEROVOL CONVERGENCE")
print("=" * 78)

settle_dt = Date(1, 12, 2019)  # CHANGED
rate = 0.05
discount_curve = FlatDiscountCurve(settle_dt, rate, FrequencyTypes.ANNUAL)

# Bond details
issue_dt = Date(1, 9, 2015)
maturity_dt = Date(1, 9, 2025)
coupon = 0.06
freq_type = FrequencyTypes.ANNUAL
dc_type = DayCountTypes.ACT_ACT_ICMA
bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

# Option Details
expiry_dt = settle_dt.add_tenor("18m")  # Date(1, 12, 2021)
#    print("EXPIRY:", expiry_dt)

df_expiry = discount_curve.df(expiry_dt)
spot_clean_value = bond.clean_price_from_discount_curve(settle_dt, discount_curve)
fwd_clean_value = bond.clean_price_from_discount_curve(expiry_dt, discount_curve)
#    print("BOND SpotCleanBondPx", spot_clean_value)
#    print("BOND FwdCleanBondPx", fwd_clean_value)
#    print("BOND Accrued:", bond.accrued_int)

spot_clean_value = bond.clean_price_from_discount_curve(settle_dt, discount_curve)

print(
    "STRIKE",
    "STEPS",
    "CALL_INT",
    "CALL_INT_PV",
    "CALL_EUR",
    "CALL_AMER",
    "PUT_INT",
    "PUT_INT_PV",
    "PUT_EUR",
    "PUT_AMER",
)

strike_prices = [90, 100, 110, 120]

for strike_price in strike_prices:

    call_intrinsic = max(spot_clean_value - strike_price, 0)
    put_intrinsic = max(strike_price - spot_clean_value, 0)
    call_intrinsic_pv = max(fwd_clean_value - strike_price, 0) * df_expiry
    put_intrinsic_pv = max(strike_price - fwd_clean_value, 0) * df_expiry

    num_time_steps = range(100, 1000, 200)

    for num_steps in num_time_steps:

        sigma = 0.0000001
        model = BDTTree(sigma, num_steps)

        opt_type = OptionTypes.EUROPEAN_CALL
        bond_option1 = BondOption(bond, expiry_dt, strike_price, opt_type)
        v1 = bond_option1.value(settle_dt, discount_curve, model)

        opt_type = OptionTypes.AMERICAN_CALL
        bond_option2 = BondOption(bond, expiry_dt, strike_price, opt_type)
        v2 = bond_option2.value(settle_dt, discount_curve, model)

        opt_type = OptionTypes.EUROPEAN_PUT
        bond_option3 = BondOption(bond, expiry_dt, strike_price, opt_type)
        v3 = bond_option3.value(settle_dt, discount_curve, model)

        opt_type = OptionTypes.AMERICAN_PUT
        bond_option4 = BondOption(bond, expiry_dt, strike_price, opt_type)
        v4 = bond_option4.value(settle_dt, discount_curve, model)

        print(
            strike_price,
            num_steps,
            call_intrinsic,
            call_intrinsic_pv,
            v1,
            v2,
            put_intrinsic,
            put_intrinsic_pv,
            v3,
            v4,
        )

# ============================================================================
# 2. BOND OPTION
# ============================================================================
# What this section demonstrates:
# Values cash flows directly from discount factors and includes accrued interest.
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("2. BOND OPTION")
print("=" * 78)

settle_dt = Date(1, 12, 2019)
issue_dt = Date(1, 12, 2018)
maturity_dt = settle_dt.add_tenor("10Y")
coupon = 0.05
freq_type = FrequencyTypes.SEMI_ANNUAL
dc_type = DayCountTypes.ACT_ACT_ICMA
bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

t_mat = (maturity_dt - settle_dt) / G_DAYS_IN_YEAR
times = np.linspace(0, t_mat, 20)
dates = settle_dt.add_years(times)
dfs = np.exp(-0.05 * times)
discount_curve = DiscountCurve(settle_dt, dates, dfs)

expiry_dt = settle_dt.add_tenor("18m")
strike_price = 105.0

strikes = [80, 100, 120]

opt_type = OptionTypes.EUROPEAN_CALL

print("LABEL", "VALUE")

price = bond.dirty_price_from_discount_curve(settle_dt, discount_curve)
print("Fixed Income Price:", price)

num_time_steps = 100

print("OPTION TYPE AND MODEL", "STRIKE", "VALUE")

for strike_price in strikes:

    sigma = 0.20

    bond_option = BondOption(bond, expiry_dt, strike_price, opt_type)
    model = BDTTree(sigma, num_time_steps)
    v = bond_option.value(settle_dt, discount_curve, model)
    print("EUROPEAN CALL - BK", strike_price, v)

for strike_price in strikes:

    sigma = 0.20

    bond_option = BondOption(bond, expiry_dt, strike_price, opt_type)
    model = BDTTree(sigma, num_time_steps)
    v = bond_option.value(settle_dt, discount_curve, model)
    print("EUROPEAN CALL - BK", strike_price, v)

opt_type = OptionTypes.AMERICAN_CALL

price = bond.dirty_price_from_discount_curve(settle_dt, discount_curve)
print("LABEL", "VALUE")
print("Fixed Income Price:", price)

print("OPTION TYPE AND MODEL", "STRIKE", "VALUE")

for strike_price in strikes:

    sigma = 0.20

    bond_option = BondOption(bond, expiry_dt, strike_price, opt_type)
    model = BDTTree(sigma, num_time_steps)
    v = bond_option.value(settle_dt, discount_curve, model)
    print("AMERICAN CALL - BK", strike_price, v)

for strike_price in strikes:

    sigma = 0.20

    bond_option = BondOption(bond, expiry_dt, strike_price, opt_type)
    model = BDTTree(sigma, num_time_steps)
    v = bond_option.value(settle_dt, discount_curve, model)
    print("AMERICAN CALL - BK", strike_price, v)

opt_type = OptionTypes.EUROPEAN_PUT

price = bond.dirty_price_from_discount_curve(settle_dt, discount_curve)

for strike_price in strikes:

    sigma = 0.01

    bond_option = BondOption(bond, expiry_dt, strike_price, opt_type)
    model = BDTTree(sigma, num_time_steps)
    v = bond_option.value(settle_dt, discount_curve, model)
    print("EUROPEAN PUT - BK", strike_price, v)

for strike_price in strikes:

    sigma = 0.20

    bond_option = BondOption(bond, expiry_dt, strike_price, opt_type)
    model = BDTTree(sigma, num_time_steps)
    v = bond_option.value(settle_dt, discount_curve, model)
    print("EUROPEAN PUT - BK", strike_price, v)

opt_type = OptionTypes.AMERICAN_PUT

price = bond.dirty_price_from_discount_curve(settle_dt, discount_curve)

for strike_price in strikes:

    sigma = 0.02

    bond_option = BondOption(bond, expiry_dt, strike_price, opt_type)
    model = BDTTree(sigma, num_time_steps)
    v = bond_option.value(settle_dt, discount_curve, model)
    print("AMERICAN PUT - BK", strike_price, v)

for strike_price in strikes:

    sigma = 0.20

    bond_option = BondOption(bond, expiry_dt, strike_price, opt_type)
    model = BDTTree(sigma, num_time_steps)
    v = bond_option.value(settle_dt, discount_curve, model)
    print("AMERICAN PUT - BK", strike_price, v)

# ============================================================================
# 3. BOND OPTION AMERICAN CONVERGENCE TWO
# ============================================================================
# What this section demonstrates:
# Values cash flows directly from discount factors and includes accrued interest.
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("3. BOND OPTION AMERICAN CONVERGENCE TWO")
print("=" * 78)

settle_dt = Date(1, 12, 2019)
discount_curve = FlatDiscountCurve(settle_dt, 0.05, FrequencyTypes.CONTINUOUS)

# Bond details
issue_dt = Date(1, 9, 2014)
maturity_dt = Date(1, 9, 2025)
coupon = 0.05
freq_type = FrequencyTypes.ANNUAL
dc_type = DayCountTypes.ACT_ACT_ICMA
bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)
expiry_dt = settle_dt.add_tenor("18m")

spot_value = bond.dirty_price_from_discount_curve(settle_dt, discount_curve)
print("LABEL", "VALUE")
print("BOND PRICE", spot_value)

print("TIME", "N", "EUR_CALL", "AMER_CALL", "EUR_PUT", "AMER_PUT")

sigma = 0.2
model = BDTTree(sigma)
k = 101.0

vec_ec = []
vec_ac = []
vec_ep = []
vec_ap = []

if True:
    k = 100.0
    bk_model = BDTTree(sigma, 100)
    euro_call_bond_option = BondOption(
        bond, expiry_dt, k, OptionTypes.EUROPEAN_CALL
    )

    v_ec = euro_call_bond_option.value(settle_dt, discount_curve, model)
    print("LABEL", "VALUE")
    print("OPTION", v_ec)

num_steps_vector = range(100, 100, 1)  # should be 100-400

for num_steps in num_steps_vector:

    bk_model = BDTTree(sigma, num_steps)

    start = time.time()

    euro_call_bond_option = BondOption(
        bond, expiry_dt, k, OptionTypes.EUROPEAN_CALL
    )
    v_ec = euro_call_bond_option.value(settle_dt, discount_curve, bk_model)

    amer_call_bond_option = BondOption(
        bond, expiry_dt, k, OptionTypes.AMERICAN_CALL
    )
    v_ac = amer_call_bond_option.value(settle_dt, discount_curve, bk_model)

    euro_put_bond_option = BondOption(bond, expiry_dt, k, OptionTypes.EUROPEAN_PUT)
    v_ep = euro_put_bond_option.value(settle_dt, discount_curve, bk_model)

    amer_put_bond_option = BondOption(bond, expiry_dt, k, OptionTypes.AMERICAN_PUT)
    v_ap = amer_put_bond_option.value(settle_dt, discount_curve, bk_model)

    end = time.time()
    period = end - start

    print(period, num_steps, v_ec, v_ac, v_ep, v_ap)

    vec_ec.append(v_ec)
    vec_ac.append(v_ac)
    vec_ep.append(v_ep)
    vec_ap.append(v_ap)

if PLOT_GRAPHS:

    plt.figure()
    plt.plot(num_steps_vector, vec_ec, label="European Call")
    plt.legend()

    plt.figure()
    plt.plot(num_steps_vector, vec_ac, label="American Call")
    plt.legend()

    plt.figure()
    plt.plot(num_steps_vector, vec_ep, label="European Put")
    plt.legend()

    plt.figure()
    plt.plot(num_steps_vector, vec_ap, label="American Put")
    plt.legend()

