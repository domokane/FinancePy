import pytest
import numpy as np
from financepy.products.rates.cross_currency_basis_swap import CrossCurrencyBasisSwap
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.date import Date
from financepy.utils.global_types import SwapTypes
from financepy.utils.error import FinError


@pytest.mark.parametrize("direction", [SwapTypes.PAY, SwapTypes.RECEIVE])
def test_cross_currency_par_swap_zero_value_at_inception(direction):
    start = Date(5, 1, 2026)
    swap = CrossCurrencyBasisSwap(
        start, "2Y", "EURUSD", direction,
        domestic_notional=1_200_000, foreign_notional=1_000_000
    )
    domestic=FlatDiscountCurve(start, .045)
    foreign=FlatDiscountCurve(start, .03)
    npv=swap.value(start, 1.2, domestic, foreign)
    assert npv == pytest.approx(0, abs=1e-6)


def test_cross_currency_foreign_spread_par_solver():
    start = Date(5, 1, 2026)
    swap = CrossCurrencyBasisSwap(
        start, "3Y", "EURUSD", SwapTypes.PAY,
        1_200_000, 1_000_000, domestic_spread=.0015,
        foreign_spread=0.0,
    )
    d=FlatDiscountCurve(start,.045)
    f=FlatDiscountCurve(start,.03)
    basis=swap.par_foreign_spread(start,1.2,d,f)
    assert 0.0 < basis < 0.01
    swap.foreign_leg.spread=basis
    assert swap.value(start,1.2,d,f)==pytest.approx(0, abs=1e-6)


def test_cross_currency_pv_matches_individual_legs_without_par_exchange():
    start=Date(5,1,2026)
    swap=CrossCurrencyBasisSwap(
        start, "1Y", "EURUSD", SwapTypes.PAY,
        1_200_000, 1_000_000,
        domestic_spread=.001, foreign_spread=.0004,
        exchange_principal=False,
    )
    d=FlatDiscountCurve(start,.045)
    f=FlatDiscountCurve(start,.03)
    got=swap.value(start,1.2,d,f)
    expected=(swap.domestic_leg.value(start,d,d)
              +1.2*swap.foreign_leg.value(start,f,f))
    assert got == pytest.approx(expected, abs=1e-9)


def test_cross_currency_initial_exchange_excluded_after_effective_date():
    start=Date(5,1,2026)
    valuation=Date(6,1,2026)
    swap=CrossCurrencyBasisSwap(
        start, "1Y", "EURUSD", SwapTypes.PAY,
        1_000_000, 900_000,
        exchange_principal=True,
    )
    d=FlatDiscountCurve(start,.04)
    f=FlatDiscountCurve(start,.02)
    expected=(swap.domestic_leg.value(valuation,d,d)
              +1.17*swap.foreign_leg.value(valuation,f,f))
    assert swap.value(valuation,1.17,d,f)==pytest.approx(expected, abs=1e-9)


@pytest.mark.parametrize("domestic_notional,foreign_notional", [(0,1), (1,0), (-1,1)])
def test_cross_currency_rejects_invalid_notionals(domestic_notional,foreign_notional):
    with pytest.raises(FinError):
        CrossCurrencyBasisSwap(Date(5,1,2026),"1Y","EURUSD",SwapTypes.PAY,
                               domestic_notional,foreign_notional)


def test_cross_currency_rejects_invalid_fx_quote():
    start=Date(5,1,2026)
    swap=CrossCurrencyBasisSwap(start,"1Y","EURUSD",SwapTypes.PAY,1_000_000,800_000)
    curve=FlatDiscountCurve(start,.04)
    with pytest.raises(FinError):
        swap.value(start,0,curve,curve)


def test_par_spread_rejects_matured_swap():
    start=Date(5,1,2026)
    swap=CrossCurrencyBasisSwap(start,"6M","EURUSD",SwapTypes.PAY,1_000_000,800_000)
    curve=FlatDiscountCurve(start,.04)
    with pytest.raises(FinError, match="No outstanding"):
        swap.par_foreign_spread(Date(8,1,2027),1.1,curve,curve)

