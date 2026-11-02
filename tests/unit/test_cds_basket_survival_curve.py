"""The basket survival curve built by the one-factor Gaussian basket pricer."""

import numpy as np
import pytest

from financepy.market.curves.cds_curve import CDSCurve
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.products.credit.cds import CDS
from financepy.products.credit.cds_basket import CDSBasket
from financepy.utils.date import Date

value_dt = Date(2, 1, 2026)
maturity_dt = Date(20, 12, 2030)
recovery = 0.4
libor_curve = FlatDiscountCurve(value_dt, 0.03)


def _issuer_curve(spread):
    return CDSCurve(value_dt, [CDS(value_dt, maturity_dt, spread)], libor_curve, recovery)


########################################################################################


@pytest.mark.parametrize("spread", [0.01, 0.05])
@pytest.mark.parametrize("beta", [0.0, 0.7])
def test_single_name_basket_reprices_the_cds(spread, beta):

    # A first-to-default basket on one credit is that credit's CDS, whatever the
    # correlation. Its breakeven spread must equal the CDS par spread.
    issuer_curve = _issuer_curve(spread)
    basket = CDSBasket(value_dt, maturity_dt)

    par_spread = basket.cds_contract.par_spread(value_dt, issuer_curve)
    basket_spread = basket.value_1f_gaussian_homo(
        value_dt, 1, [issuer_curve], np.array([beta]), libor_curve, 50
    )[3]

    assert basket_spread == pytest.approx(par_spread, rel=1e-5)


def test_independent_first_to_default_matches_product_of_survival_curves():

    # With zero correlation the first-to-default survival probability is the
    # product of the single-name survival probabilities.
    issuer_curves = [_issuer_curve(0.01), _issuer_curve(0.02)]
    basket = CDSBasket(value_dt, maturity_dt)

    product_curve = CDSCurve(value_dt, [], libor_curve, recovery)
    times = np.linspace(0.0, 5.0, 2001)
    product_curve._times = times
    product_curve._qs = np.array(
        [np.prod([curve.survival_prob(t) for curve in issuer_curves]) for t in times]
    )
    reference_spread = basket.cds_contract.par_spread(value_dt, product_curve)

    basket_spread = basket.value_1f_gaussian_homo(
        value_dt, 1, issuer_curves, np.zeros(2), libor_curve, 50
    )[3]

    assert basket_spread == pytest.approx(reference_spread, rel=1e-5)
