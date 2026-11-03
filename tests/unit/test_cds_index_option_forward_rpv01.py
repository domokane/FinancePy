"""The forward RPV01 inside the Anderson CDS index option pricer."""

from math import exp

import pytest

from financepy.market.curves.cds_curve import CDSCurve
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.products.credit.cds import CDS
from financepy.products.credit.cds_index_option import CDSIndexOption
from financepy.utils.date import Date
from financepy.utils.day_count import DayCount
from financepy.utils.helpers import times_from_dates

value_dt = Date(2, 1, 2026)
maturity_dt = Date(20, 12, 2030)
index_cpn = 0.01
recovery = 0.4
libor_curve = FlatDiscountCurve(value_dt, 0.03)


def _forward_rpv01_from_pricer(option, spread):
    # With zero volatility the Anderson integrand is (s - c) * fwd_rpv01(s) at s = x
    # for every quadrature node, so the expected value returns the forward RPV01.
    h, _ = option._calc_index_payer_option_price(
        value_dt, spread, 0.0, index_cpn, 0.0, libor_curve, recovery
    )
    return h / (spread - index_cpn)


def _forward_rpv01_by_hand(option, spread):
    # Every coupon of the forward CDS, survival weighted with a flat hazard rate
    # and discounted from expiry, less the premium accrued from the previous
    # coupon date to expiry, which the protection buyer owes at the start.
    contract = option.cds_contract
    expiry_dt = option.expiry_dt
    dc_type = libor_curve.curve_dc_type
    lgd = 1.0 - recovery
    df_to_expiry = libor_curve.df(expiry_dt)
    total = 0.0
    for accrual_factor, payment_dt in zip(contract.accrual_factors, contract.payment_dts):
        t = times_from_dates(expiry_dt, payment_dt, dc_type)
        total += accrual_factor * exp(-spread * t / lgd) * libor_curve.df(payment_dt) / df_to_expiry
    day_count = DayCount(option.accrual_dc_type)
    accrued = day_count.year_frac(contract.accrual_start_dts[0], expiry_dt)[0]
    return total - max(accrued, 0.0)


########################################################################################


@pytest.mark.parametrize("expiry_dt", [Date(20, 12, 2026), Date(15, 2, 2027), Date(19, 6, 2027)])
@pytest.mark.parametrize("spread", [0.005, 0.015, 0.04])
def test_forward_rpv01_counts_every_coupon_and_the_accrued_premium(expiry_dt, spread):

    option = CDSIndexOption(expiry_dt, maturity_dt, index_cpn, 0.012)
    contract = option.cds_contract

    # The forward contract starts at expiry, so its first payment is after expiry
    assert contract.payment_dts[0] > expiry_dt

    from_pricer = _forward_rpv01_from_pricer(option, spread)
    by_hand = _forward_rpv01_by_hand(option, spread)

    assert from_pricer == pytest.approx(by_hand, rel=1e-6)


def test_accrued_premium_is_not_negative_when_expiry_is_a_coupon_date():

    # 20 December 2026 is a Sunday, so the first accrual start is adjusted to the
    # Monday after expiry. Nothing has accrued, and the pricer must not add a
    # negative accrual to the forward RPV01.
    expiry_dt = Date(20, 12, 2026)
    option = CDSIndexOption(expiry_dt, maturity_dt, index_cpn, 0.012)
    contract = option.cds_contract
    assert contract.accrual_start_dts[0] > expiry_dt

    spread = 0.015
    from_pricer = _forward_rpv01_from_pricer(option, spread)
    by_hand = _forward_rpv01_by_hand(option, spread)
    assert from_pricer == pytest.approx(by_hand, rel=1e-6)


def test_anderson_close_to_adjusted_black():

    expiry_dt = Date(20, 12, 2026)
    issuer_curves = [
        CDSCurve(value_dt, [CDS(value_dt, maturity_dt, 0.012)], libor_curve, recovery)
        for _ in range(10)
    ]
    index_curve = issuer_curves[0]
    sigma = 0.3

    for strike in [0.008, 0.012, 0.016]:
        option = CDSIndexOption(expiry_dt, maturity_dt, index_cpn, strike)
        v_pay_anderson = option.value_anderson(value_dt, issuer_curves, recovery, sigma)[0]
        v_pay_black = option.value_adjusted_black(value_dt, index_curve, recovery, libor_curve, sigma)[0]
        assert v_pay_anderson == pytest.approx(v_pay_black, rel=0.02)
