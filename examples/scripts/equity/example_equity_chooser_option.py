# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.

import numpy as np
from financepy.products.equity.equity_chooser_option import EquityChooserOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.date import Date

# ============================================================================
# FINANCEPY EXAMPLES - EquityChooserOption
# ============================================================================



########################################################################################




########################################################################################




########################################################################################




########################################################################################



# %% chooser_value_as_function_of_stock_price
value_dt = Date(1, 1, 2027)
stock_price = 100.0
volatility = 0.20
interest_rate = 0.04
dividend_yield = 0.02
discount_curve = FlatDiscountCurve(value_dt, interest_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

choose_dt = Date(1, 6, 2027)
call_K = 100.0
put_K = 100.0
call_exp_dt = Date(1, 1, 2028)
put_exp_dt = Date(1, 1, 2028)
model = BlackScholes(volatility)

chooser = EquityChooserOption(choose_dt, call_exp_dt, put_exp_dt, call_K, put_K)

s_vec = np.linspace(50, 150, 100)
v_vec = chooser.value(value_dt, s_vec, discount_curve, dividend_curve, model)

# ============================================================================
# 1. EQUITY CHOOSER OPTION HAUG
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.

print("\n" + "=" * 78)
print("1. EQUITY CHOOSER OPTION HAUG")
print("=" * 78)

"""Following example in Haug Page 130"""

value_dt = Date(1, 1, 2015)
choose_dt = Date(2, 4, 2015)
call_expiry_dt = Date(1, 7, 2015)
put_expiry_dt = Date(2, 8, 2015)
call_strike = 55.0
put_strike = 48.0
stock_price = 50.0
volatility = 0.35
interest_rate = 0.10
dividend_yield = 0.05

model = BlackScholes(volatility)
discount_curve = FlatDiscountCurve(value_dt, interest_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

chooser_option = EquityChooserOption(choose_dt, call_expiry_dt, put_expiry_dt, call_strike, put_strike)

v = chooser_option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

v_mc = chooser_option.value_mc(value_dt, stock_price, discount_curve, dividend_curve, model, 20000)

v_haug = 6.0508
print("", "", "", "", "", "")
print("FINANCEPY", v, "HAUG", v_haug, "MC", v_mc)

# ============================================================================
# 2. EQUITY CHOOSER OPTION MATLAB
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.

print("\n" + "=" * 78)
print("2. EQUITY CHOOSER OPTION MATLAB")
print("=" * 78)

"""https://fr.mathworks.com/help/fininst/chooserbybls.html"""

value_dt = Date(1, 6, 2007)
choose_date = Date(31, 8, 2007)
call_expiry_dt = Date(2, 12, 2007)
put_expiry_dt = Date(2, 12, 2007)
call_strike = 60.0
put_strike = 60.0
stock_price = 50.0
volatility = 0.20
interest_rate = 0.10
dividend_yield = 0.05

model = BlackScholes(volatility)

discount_curve = FlatDiscountCurve(value_dt, interest_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

chooser_option = EquityChooserOption(choose_date, call_expiry_dt, put_expiry_dt, call_strike, put_strike)

v = chooser_option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

v_mc = chooser_option.value_mc(value_dt, stock_price, discount_curve, dividend_curve, model, 20000)

v_matlab = 8.9308
print("", "", "", "", "", "")
print("FINANCEPY", v, "MATLAB", v_matlab, "MC", v_mc)

# ============================================================================
# 3. EQUITY CHOOSER OPTION DERIVICOM
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.

print("\n" + "=" * 78)
print("3. EQUITY CHOOSER OPTION DERIVICOM")
print("=" * 78)

"""http://derivicom.com/support/finoptionsxl/index.html?complex_chooser.htm"""

value_dt = Date(1, 1, 2007)
choose_date = Date(1, 2, 2007)
call_expiry_dt = Date(1, 4, 2007)
put_expiry_dt = Date(1, 5, 2007)
call_strike = 40.0
put_strike = 35.0
stock_price = 38.0
volatility = 0.20
interest_rate = 0.08
dividend_yield = 0.0625

model = BlackScholes(volatility)
discount_curve = FlatDiscountCurve(value_dt, interest_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

chooser_option = EquityChooserOption(choose_date, call_expiry_dt, put_expiry_dt, call_strike, put_strike)

v = chooser_option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

v_mc = chooser_option.value_mc(value_dt, stock_price, discount_curve, dividend_curve, model, 20000)

v_derivicom = 1.0989
print("", "", "", "", "", "")
print("FINANCEPY", v, "DERIVICOM", v_derivicom, "MC", v_mc)

