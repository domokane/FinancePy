"""Basket PV01 includes the real first scheduled premium accrual."""

import numpy as np
import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.market.curves.cds_curve import CDSCurve
from financepy.products.credit.cds import CDS
from financepy.products.credit.cds_basket import CDSBasket
from financepy.utils.date import Date
from financepy.utils.day_count import DayCount, DayCountTypes


@pytest.mark.parametrize("basis", [DayCountTypes.ACT_360, DayCountTypes.ACT_365F])
def test_no_default_pv01_includes_first_coupon(basis):
    """No-default paths equal the discounted sum of all scheduled accruals."""
    today = Date(2, 1, 2026)
    maturity = Date(20, 12, 2026)
    curve = FlatDiscountCurve(today, .03)
    basket = CDSBasket(today, maturity, accrual_dc_type=basis)
    issuer = CDSCurve(today, [CDS(today, maturity, .01)], curve, 0.4)
    pv01, protection = basket.value_legs_mc(today, 1, np.array([[100.0, 100.0]]), [issuer], curve)
    contract = basket.cds_contract
    starts = [contract.accrual_start_dts[0]] + contract.payment_dts[:-1]
    day_count = DayCount(basis)
    expected = sum(day_count.year_frac(start, end)[0] * curve.df(end) for start, end in zip(starts, contract.payment_dts))
    assert pv01 == pytest.approx(expected)
    assert protection == 0.0
    assert expected > sum(day_count.year_frac(start, end)[0] * curve.df(end) for start, end in zip(contract.payment_dts[:-1], contract.payment_dts[1:]))


def test_first_period_default_accrues_from_actual_first_start():
    """A first-period default owes the accrued premium, not zero."""
    today = Date(2, 1, 2026)
    maturity = Date(20, 12, 2026)
    curve = FlatDiscountCurve(today, .03)
    basket = CDSBasket(today, maturity, accrual_dc_type=DayCountTypes.ACT_365F)
    issuer = CDSCurve(today, [CDS(today, maturity, .01)], curve, 0.4)
    tau = 0.02
    assert tau < (basket.cds_contract.payment_dts[0] - today) / 365.0
    pv01, protection = basket.value_legs_mc(today, 1, np.array([[tau]]), [issuer], curve)
    start = (basket.cds_contract.accrual_start_dts[0] - today) / 365.0
    assert pv01 == pytest.approx((tau - start) * curve.df_t(tau))
    assert protection == pytest.approx((1.0 - issuer.recovery_rate) * curve.df_t(tau))
