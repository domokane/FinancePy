# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.

from financepy.products.equity.equity_cliquet_option import EquityCliquetOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes

# ============================================================================
# FINANCEPY EXAMPLES - EquityCliquetOption
# ============================================================================



########################################################################################




########################################################################################

# ============================================================================
# 1. EQUITY CLIQUET OPTION
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.

print("\n" + "=" * 78)
print("1. EQUITY CLIQUET OPTION")
print("=" * 78)

start_dt = Date(1, 1, 2014)
final_expiry_dt = Date(1, 1, 2017)
freq_type = FrequencyTypes.QUARTERLY
opt_type = OptionTypes.EUROPEAN_CALL

cliquet_option = EquityCliquetOption(start_dt, final_expiry_dt, opt_type, freq_type)

value_dt = Date(1, 1, 2015)
stock_price = 100.0
volatility = 0.20
interest_rate = 0.05
dividend_yield = 0.02
model = BlackScholes(volatility)
discount_curve = FlatDiscountCurve(value_dt, interest_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

v = cliquet_option.value(
    value_dt, stock_price, discount_curve, dividend_curve, model
)

print("LABEL", "VALUE")
print("FINANCEPY", v)

