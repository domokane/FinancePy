# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.

from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes


from financepy.utils.global_types import OptionTypes
from financepy.utils.global_types import BlackScholesTypes

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.products.equity.equity_american_option import EquityAmericanOption


PLOT_GRAPHS = False



import matplotlib.pyplot as plt

# ============================================================================
# FINANCEPY EXAMPLES - BlackScholes
# ============================================================================

# TODO Complete output of results to log files

########################################################################################




########################################################################################

# ============================================================================
# 1. BARONE EDESI
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("1. BARONE EDESI")
print("=" * 78)

value_dt = Date(8, 5, 2015)
expiry_dt = Date(15, 1, 2016)

strike_price = 130.0
stock_price = 127.62
volatility = 0.20
interest_rate = 0.001
dividend_yield = 0.0163

opt_type = OptionTypes.AMERICAN_CALL
eu_option_type = OptionTypes.EUROPEAN_CALL

am_option = EquityAmericanOption(expiry_dt, strike_price, opt_type)

ameu_option = EquityAmericanOption(expiry_dt, strike_price, eu_option_type)

eu_option = EquityVanillaOption(expiry_dt, strike_price, eu_option_type)

discount_curve = FlatDiscountCurve(
    value_dt, interest_rate, FrequencyTypes.CONTINUOUS, DayCountTypes.ACT_365F
)

dividend_curve = FlatDiscountCurve(
    value_dt, dividend_yield, FrequencyTypes.CONTINUOUS, DayCountTypes.ACT_365F
)

num_steps_per_year = 400

model_tree = BlackScholes(volatility, BlackScholesTypes.CRR_TREE, num_steps_per_year)

v = am_option.value(
     value_dt, stock_price, discount_curve, dividend_curve, model_tree
 )
assert round(v, 4) == 6.8398

model_approx = BlackScholes(volatility, BlackScholesTypes.BARONE_ADESI)

v = am_option.value(
     value_dt, stock_price, discount_curve, dividend_curve, model_approx
)

assert round(v, 4) == 6.8277

v = ameu_option.value(
     value_dt, stock_price, discount_curve, dividend_curve, model_tree
)

assert round(v, 4) == 6.7512

v = eu_option.value(
     value_dt, stock_price, discount_curve, dividend_curve, model_tree
)

assert round(v, 4) == 6.7493

# ============================================================================
# 2. BLACK SCHOLES
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("2. BLACK SCHOLES")
print("=" * 78)

value_dt = Date(8, 5, 2015)
expiry_dt = Date(15, 1, 2016)

strike_price = 130.0
stock_price = 127.62
volatility = 0.20
interest_rate = 0.001
dividend_yield = 0.0163

opt_type = OptionTypes.AMERICAN_CALL
eu_option_type = OptionTypes.EUROPEAN_CALL

# Pure American Option TREE
am_option = EquityAmericanOption(expiry_dt, strike_price, opt_type)

# American with European style exercise TREE
ameu_option = EquityAmericanOption(expiry_dt, strike_price, eu_option_type)

# European Option and European Exercise so Black Scholes
eu_option = EquityVanillaOption(expiry_dt, strike_price, eu_option_type)

discount_curve = FlatDiscountCurve(
    value_dt,
    interest_rate,
    FrequencyTypes.CONTINUOUS,
)

dividend_curve = FlatDiscountCurve(
    value_dt,
    dividend_yield,
    FrequencyTypes.CONTINUOUS,
)

am_tree_value = []
am_baw_value = []
eu_tree_value = []
eu_anal_value = []
volatility = 0.20

num_steps_per_year = range(5, 200, 1)

print(
    "STEPS PER YEAR",
    "AMERICAN_TREE",
    "AMERICAN_BAW",
    "EUROPEAN_TREE",
    "EUROPEAN_BS",
)

for num_steps in num_steps_per_year:

    model_tree = BlackScholes(volatility, BlackScholesTypes.CRR_TREE, num_steps)
    model_anal = BlackScholes(volatility, BlackScholesTypes.ANALYTICAL)
    model_BAW = BlackScholes(volatility, BlackScholesTypes.BARONE_ADESI)

    v_am = am_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model_tree
    )

    v_eu = ameu_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model_tree
    )

    v_bs = eu_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model_anal
    )

    v_am_baw = am_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model_BAW
    )

    am_tree_value.append(v_am)
    eu_tree_value.append(v_eu)
    eu_anal_value.append(v_bs)
    am_baw_value.append(v_am_baw)

    print(num_steps, v_am, v_am_baw, v_eu, v_bs)

if PLOT_GRAPHS:
    plt.title("American Option Price Convergence Analysis")
    plt.plot(num_steps_per_year, am_tree_value, label="American Tree")
    plt.plot(num_steps_per_year, am_baw_value, label="American BAW")
    plt.plot(num_steps_per_year, eu_tree_value, label="European Tree")
    plt.plot(num_steps_per_year, eu_anal_value, label="European Anal", lw=2)
    plt.xlabel("Num Steps")
    plt.ylabel("Value")
    plt.legend()
    plt.show()

