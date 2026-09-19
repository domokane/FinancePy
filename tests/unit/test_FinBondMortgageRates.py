"""Mortgage payments and amortization near the zero-rate limit."""

from decimal import Decimal, localcontext
import math

import pytest

from financepy.products.bonds.bond_mortgage import BondMortgage, BondMortgageTypes
from financepy.utils.date import Date
from financepy.utils.error import FinError
from financepy.utils.frequency import FrequencyTypes, annual_frequency


@pytest.mark.parametrize("rate", [0.0, 1e-15, -1e-15, 1e-12, -1e-12, 0.035, -0.01])
@pytest.mark.parametrize("years", [1, 30])
@pytest.mark.parametrize("frequency", [FrequencyTypes.ANNUAL,
                                      FrequencyTypes.SEMI_ANNUAL,
                                      FrequencyTypes.QUARTERLY,
                                      FrequencyTypes.MONTHLY])
def test_repayment_matches_discounted_cash_flows(rate, years, frequency):
    """A constant payment has PV equal to principal and amortizes the debt."""
    principal = 120000.0
    mortgage = BondMortgage(Date(1, 1, 2025), Date(1, 1, 2025 + years),
                            principal, freq_type=frequency)
    periods = len(mortgage.schedule.adjusted_dts) - 1
    with localcontext() as context:
        context.prec = 80
        periodic_rate = Decimal.from_float(rate) / Decimal(
            int(annual_frequency(frequency)))
        discount = Decimal(1) / (1 + periodic_rate)
        weight = Decimal(1)
        present_value = Decimal(0)
        for _ in range(periods):
            weight *= discount
            present_value += weight
        expected = float(Decimal.from_float(principal) / present_value)

    assert mortgage.repayment_amount(rate) == pytest.approx(expected, rel=2e-12)
    mortgage.generate_flows(rate, BondMortgageTypes.REPAYMENT)
    assert len(mortgage.total_flows) == periods + 1
    assert mortgage.total_flows[1:] == pytest.approx([expected] * periods, rel=2e-12)
    assert mortgage.principal_remaining[-1] == pytest.approx(0, abs=principal * 2e-10)
    assert math.fsum(mortgage.principal_flows) == pytest.approx(principal, rel=2e-10)
    if rate == 0.0:
        assert mortgage.interest_flows == [0] * (periods + 1)
        assert mortgage.principal_flows[1:] == pytest.approx(
            [principal / periods] * periods)


@pytest.mark.parametrize("rate", [0.0, 1e-12, -1e-12, 0.035])
def test_interest_only_flows_are_unchanged(rate):
    """Interest-only cash flows do not use the repayment annuity formula."""
    principal = 120000.0
    mortgage = BondMortgage(Date(1, 1, 2025), Date(1, 1, 2055), principal)
    mortgage.generate_flows(rate, BondMortgageTypes.INTEREST_ONLY)
    assert mortgage.principal_remaining == [principal] * 361
    assert mortgage.principal_flows == [0.0] * 361
    assert mortgage.total_flows[1:] == pytest.approx([principal * rate / 12] * 360)


@pytest.mark.parametrize("rate", [-12.0, -13.0])
def test_invalid_periodic_rate(rate):
    """Positive discount factors require a periodic rate greater than -1."""
    mortgage = BondMortgage(Date(1, 1, 2025), Date(1, 1, 2055), 120000.0)
    with pytest.raises(FinError, match="Periodic mortgage rate"):
        mortgage.repayment_amount(rate)
