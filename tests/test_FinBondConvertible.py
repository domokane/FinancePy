###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

import numpy as np

from financepy.products.bonds.bond_convertible import BondConvertible
from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat

maturity_date = Date(15, 3, 2022)
coupon = 0.0575
freq_type = FrequencyTypes.SEMI_ANNUAL
start_convert_date = Date(31, 12, 2003)
conversion_ratio = 3.84615  # adjust for face

call_dates = [Date(20, 3, 2007),
              Date(15, 3, 2012),
              Date(15, 3, 2017)]
call_price = 1100
call_prices = np.array([call_price, call_price, call_price])

put_dates = [Date(20, 3, 2007),
             Date(15, 3, 2012),
             Date(15, 3, 2017)]

putPrice = 90.0
put_prices = np.array([putPrice, putPrice, putPrice])
accrualBasis = DayCountTypes.ACT_365F

bond = BondConvertible(maturity_date,
                       coupon,
                       freq_type,
                       start_convert_date,
                       conversion_ratio,
                       call_dates,
                       call_prices,
                       put_dates,
                       put_prices,
                       accrualBasis)

settlement_date = Date(31, 12, 2003)
stock_price = 28.5
stock_volatility = 0.370
dividend_dates = [Date(20, 3, 2007),
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
                  Date(15, 3, 2022)]
rate = 0.04

discount_curve = DiscountCurveFlat(settlement_date,
                                   rate,
                                   FrequencyTypes.CONTINUOUS)
credit_spread = 0.00
recovery_rate = 0.40


def test_calls_or_puts():

    dividend_yields = [0.00] * 16

    num_steps_per_year = 5

    res = bond.value(settlement_date,
                     stock_price,
                     stock_volatility,
                     dividend_dates,
                     dividend_yields,
                     discount_curve,
                     credit_spread,
                     recovery_rate,
                     num_steps_per_year)

    assert round(res['cbprice'], 4) == 197.3373
    assert round(res['bond'], 4) == 123.5343
    assert round(res['delta'], 4) == 3.1938
    assert round(res['gamma'], 4) == 0.0153
    assert round(res['theta'], 4) == 39.6423

    num_steps_per_year = 20
    res = bond.value(settlement_date,
                     stock_price,
                     stock_volatility,
                     dividend_dates,
                     dividend_yields,
                     discount_curve,
                     credit_spread,
                     recovery_rate,
                     num_steps_per_year)

    assert round(res['cbprice'], 4) == 201.7815
    assert round(res['bond'], 4) == 123.5343
    assert round(res['delta'], 4) == 3.3634
    assert round(res['gamma'], 4) == 0.0832
    assert round(res['theta'], 4) == 53.9851


def test_dividends():
    dividend_yields = [0.02] * 16

    num_steps_per_year = 5
    res = bond.value(settlement_date,
                     stock_price,
                     stock_volatility,
                     dividend_dates,
                     dividend_yields,
                     discount_curve,
                     credit_spread,
                     recovery_rate,
                     num_steps_per_year)

    assert round(res['cbprice'], 4) == 192.133
    assert round(res['bond'], 4) == 123.5343
    assert round(res['delta'], 4) == 3.0904
    assert round(res['gamma'], 4) == 0.0157
    assert round(res['theta'], 4) == 37.8822

    num_steps_per_year = 20
    res = bond.value(settlement_date,
                     stock_price,
                     stock_volatility,
                     dividend_dates,
                     dividend_yields,
                     discount_curve,
                     credit_spread,
                     recovery_rate,
                     num_steps_per_year)

    assert round(res['cbprice'], 4) == 181.8036
    assert round(res['bond'], 4) == 123.5343
    assert round(res['delta'], 4) == 2.7362
    assert round(res['gamma'], 4) == 0.0784
    assert round(res['theta'], 4) == 44.6321
