# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path

_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause

_install_double_click_pause()
import numpy as np

import add_fp_to_path

from financepy.models.cev import CEV
from financepy.models.dupire import Dupire

########################################################################################


def test_dupire_recovers_cev():

    stock_price = 100.0

    interest_rate = 0.05
    dividend_yield = 0.02

    spot_vol = 0.20

    # Use a reasonably dense grid because Dupire requires first and
    # second derivatives of the call-price surface.
    strikes = np.linspace(
        50.0,
        150.0,
        101,
    )

    expiries = np.linspace(
        0.10,
        5.00,
        50,
    )

    # Test points deliberately lie away from the interpolation boundaries.
    test_strikes = np.array(
        [
            70.0,
            85.0,
            100.0,
            115.0,
            130.0,
        ]
    )

    test_expiries = np.array(
        [
            0.25,
            0.50,
            1.00,
            2.00,
            4.00,
        ]
    )

    betas = [
        1.0 / 3.0,
        2.0 / 3.0,
        1.0,
    ]

    print(
        "BETA",
        "T",
        "K",
        "EXPECTED",
        "DUPIRE",
        "ERROR",
    )

    for beta in betas:

        # Keep the initial local volatility at S0 equal to 20%.
        #
        # sigma_loc(S) = sigma * S^(beta - 1)
        #
        # therefore
        #
        # sigma = spot_vol * S0^(1 - beta)

        sigma = spot_vol * stock_price ** (1.0 - beta)

        cev_model = CEV(
            sigma=sigma,
            beta=beta,
        )

        # Build the CEV call-price surface.
        call_prices = np.empty(
            (
                len(expiries),
                len(strikes),
            ),
            dtype=float,
        )

        for i, t_exp in enumerate(expiries):

            for j, strike in enumerate(strikes):

                call_prices[i, j] = cev_model.call_value(
                    stock_price,
                    t_exp,
                    strike,
                    interest_rate,
                    dividend_yield,
                )

        # Recover local volatility from the option-price surface.
        dupire_model = Dupire(
            expiries,
            strikes,
            interest_rate,
            dividend_yield,
            call_prices,
        )

        for t_exp in test_expiries:

            for strike in test_strikes:

                dupire_vol = dupire_model.local_volatility(
                    strike,
                    t_exp,
                )

                # CEV percentage local volatility:
                #
                # sigma_loc(S) = sigma * S^(beta - 1)
                #
                # Dupire evaluates this at S = K.
                expected_vol = sigma * strike ** (beta - 1.0)

                error = dupire_vol - expected_vol

                assert np.isfinite(dupire_vol)
                assert abs(error) < 1.0e-2

                print(
                    beta,
                    t_exp,
                    strike,
                    expected_vol,
                    dupire_vol,
                    error,
                )


########################################################################################


test_dupire_recovers_cev()
