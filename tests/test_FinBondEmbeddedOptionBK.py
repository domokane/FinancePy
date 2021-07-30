###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes

from financepy.products.rates.ibor_swap import IborSwap

from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.products.bonds.bond import Bond
from financepy.products.bonds.bond_callable import BondEmbeddedOption
from financepy.utils.global_types import SwapTypes

from financepy.models.bk_tree import BKTree


def test_matlab_example():
    # https://fr.mathworks.com/help/fininst/optembndbybk.html
    issue_date = Date(1, 1, 2005)
    maturity_date = Date(1, 1, 2010)
    coupon = 0.0525
    freq_type = FrequencyTypes.ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_date, maturity_date, coupon, freq_type, accrual_type)
    
    valuation_date = Date(1, 1, 2007)
    settlement_date = valuation_date

    fixed_leg_type = SwapTypes.PAY
    dcType = DayCountTypes.THIRTY_E_360
    fixedFreq = FrequencyTypes.ANNUAL
    swap1 = IborSwap(settlement_date, "1Y", fixed_leg_type, 0.0350, fixedFreq, dcType)
    swap2 = IborSwap(settlement_date, "2Y", fixed_leg_type, 0.0400, fixedFreq, dcType)
    swap3 = IborSwap(settlement_date, "3Y", fixed_leg_type, 0.0450, fixedFreq, dcType)
    swaps = [swap1, swap2, swap3]
    discount_curve = IborSingleCurve(valuation_date, [], [], swaps)
    
    v = bond.clean_price_from_discount_curve(settlement_date, discount_curve)

    assert round(v, 4) == 102.0746


    call_dates = []
    call_prices = []
    put_dates = []
    put_prices = []
    putDate = Date(1, 1, 2008)
    for _ in range(0, 24):
        put_dates.append(putDate)
        put_prices.append(100)
        putDate = putDate.add_months(1)
    puttableBond = BondEmbeddedOption(issue_date, maturity_date, coupon,
                                      freq_type, accrual_type,
                                      call_dates, call_prices,
                                      put_dates, put_prices)

    sigma = 0.01  # This volatility is very small for a BK process
    a = 0.1
    num_time_steps = 100
    model = BKTree(sigma, a, num_time_steps)

    v = puttableBond.value(settlement_date, discount_curve, model)

    assert round(v['bondwithoption'], 4) == 102.3625
    assert round(v['bondpure'], 4) == 102.0746


def test_quantlib_example():
    # http://gouthamanbalaraman.com/blog/callable-bond-quantlib-python.html
    issue_date = Date(15, 9, 2010)
    maturity_date = Date(15, 9, 2022)
    coupon = 0.025
    freq_type = FrequencyTypes.QUARTERLY
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_date, maturity_date, coupon, freq_type, accrual_type)

    valuation_date = Date(16, 8, 2016)
    settlement_date = valuation_date.add_weekdays(3)

    discount_curve = DiscountCurveFlat(valuation_date, 0.035,
                                       FrequencyTypes.SEMI_ANNUAL)

    v = bond.clean_price_from_discount_curve(settlement_date, discount_curve)

    assert round(v, 4) == 94.6318


    nextCallDate = Date(15, 9, 2016)
    call_dates = [nextCallDate]
    call_prices = [100.0]
    put_dates = []
    put_prices = []
    for _ in range(1, 24):
        nextCallDate = nextCallDate.add_months(3)
        call_dates.append(nextCallDate)
        call_prices.append(100.0)
    puttableBond = BondEmbeddedOption(issue_date, maturity_date, coupon,
                                      freq_type, accrual_type,
                                      call_dates, call_prices,
                                      put_dates, put_prices)
    
    sigma = 0.12/0.035  # basis point volatility
    a = 0.03
    num_time_steps = 100
    model = BKTree(sigma, a, num_time_steps)

    v = puttableBond.value(settlement_date, discount_curve, model)

    assert round(v['bondwithoption'], 4) == 89.7614
    assert round(v['bondpure'], 4) == 95.0619
