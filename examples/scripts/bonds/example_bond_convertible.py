# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import time
import numpy as np

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date
from financepy.products.bonds.bond_convertible import BondConvertible

# ============================================================================
# FINANCEPY EXAMPLES - BondConvertible
# ============================================================================



########################################################################################




########################################################################################

# ============================================================================
# 1. BOND CONVERTIBLE
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. BOND CONVERTIBLE")
print("=" * 78)

settle_dt = Date(31, 12, 2003)
start_convert_date = Date(31, 12, 2003)
maturity_dt = Date(15, 3, 2022)
conversion_ratio = 38.4615  # adjust for face
coupon = 0.0575
freq_type = FrequencyTypes.SEMI_ANNUAL
accrual_basis = DayCountTypes.ACT_365F

call_price = 1100
call_dts = [Date(20, 3, 2007), Date(15, 3, 2012), Date(15, 3, 2017)]
call_prices = np.array([call_price, call_price, call_price])

put_price = 90
put_dts = [Date(20, 3, 2007), Date(15, 3, 2012), Date(15, 3, 2017)]
put_prices = np.array([put_price, put_price, put_price])

bond = BondConvertible(
    maturity_dt,
    coupon,
    freq_type,
    start_convert_date,
    conversion_ratio,
    call_dts,
    call_prices,
    put_dts,
    put_prices,
    accrual_basis,
)
#    print(bond)

dividend_dts = [
    Date(20, 3, 2007),
    Date(15, 3, 2008),
    Date(15, 3, 2009),
    Date(15, 3, 2010),
    Date(15, 3, 2011),
    Date(15, 3, 2012),
    Date(15, 3, 2013),
    Date(15, 3, 2014),
    Date(15, 3, 2015),
    Date(15, 3, 2016),
    Date(15, 3, 2017),
    Date(15, 3, 2018),
    Date(15, 3, 2019),
    Date(15, 3, 2020),
    Date(15, 3, 2021),
    Date(15, 3, 2022),
]

dividend_yields = [0.00] * 16
stock_price = 28.5
stock_volatility = 0.370
rate = 0.04
discount_curve = FlatDiscountCurve(settle_dt, rate, FrequencyTypes.CONTINUOUS)
credit_spread = 0.00
recovery_rate = 0.40
num_steps_per_year = 20

print("LABEL")
print("NO CALLS OR PUTS")

print("TIME", "NUMSTEPS", "PRICE")

for num_steps_per_year in [5, 10, 20]:

    start = time.time()

    res = bond.value(
        settle_dt,
        stock_price,
        stock_volatility,
        dividend_dts,
        dividend_yields,
        discount_curve,
        credit_spread,
        recovery_rate,
        num_steps_per_year,
    )

    end = time.time()
    period = end - start
    print(period, num_steps_per_year, res)

dividend_yields = [0.02] * 16
print("LABEL")
print("DIVIDENDS")

print("TIME", "NUMSTEPS", "PRICE")
for num_steps_per_year in [5, 20, 80]:
    start = time.time()
    res = bond.value(
        settle_dt,
        stock_price,
        stock_volatility,
        dividend_dts,
        dividend_yields,
        discount_curve,
        credit_spread,
        recovery_rate,
        num_steps_per_year,
    )

    end = time.time()
    period = end - start
    print(period, num_steps_per_year, res)

