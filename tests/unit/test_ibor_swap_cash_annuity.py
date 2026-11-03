# Copyright (C) 2026 Xamit Kadirbekov
# SPDX-License-Identifier: GPL-3.0-or-later

import math
import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black import Black
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.rates.ibor_swaption import IborSwaption
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.global_types import SwapTypes


FREQUENCIES = [(1, FrequencyTypes.ANNUAL), (2, FrequencyTypes.SEMI_ANNUAL),
               (4, FrequencyTypes.QUARTERLY), (12, FrequencyTypes.MONTHLY)]
VALUE_DATE = Date(15, 1, 2026)
START_DATE = Date(15, 1, 2027)


def make_swap(m, frequency, periods):
    return IborSwap(START_DATE, START_DATE.add_months(periods * (12 // m)),
                    SwapTypes.PAY, 0.04, frequency, DayCountTypes.THIRTY_E_360,
                    bd_type=BusDayAdjustTypes.NONE)


@pytest.mark.parametrize('m,frequency', FREQUENCIES)
@pytest.mark.parametrize('periods', [1, 5])
def test_zero_rate_annuity_counts_every_future_payment(m, frequency, periods):
    swap = make_swap(m, frequency, periods)
    assert len(swap.fixed_leg.payment_dts) == periods
    assert swap.cash_settled_pv01(VALUE_DATE, 0.0, frequency) == pytest.approx(periods / m)


@pytest.mark.parametrize('m,frequency', FREQUENCIES)
def test_three_payments_at_positive_rate(m, frequency):
    swap = make_swap(m, frequency, 3)
    expected = sum((1 / m) / (1 + 0.04 / m)**j for j in [1, 2, 3])
    assert swap.cash_settled_pv01(START_DATE, 0.04, frequency) == pytest.approx(expected)


@pytest.mark.parametrize('side', [SwapTypes.PAY, SwapTypes.RECEIVE])
def test_one_payment_cash_swaption_has_positive_black_value(side):
    option = IborSwaption(VALUE_DATE, START_DATE, START_DATE.add_years(1), side,
                         0.04, FrequencyTypes.ANNUAL, DayCountTypes.THIRTY_E_360,
                         notional=1e6, bd_type=BusDayAdjustTypes.NONE)
    curve = FlatDiscountCurve(VALUE_DATE, 0.03)
    # ATM Black price: F * (2*Phi(sigma*sqrt(T)/2)-1), for T=1.
    black_price = 0.04 * math.erf(0.25 / (2 * math.sqrt(2)))
    expected = 1e6 * math.exp(-0.03) * black_price / 1.04
    actual = option.cash_settled_value(VALUE_DATE, curve, 0.04, Black(0.25))
    assert actual == pytest.approx(expected, abs=0.02)
