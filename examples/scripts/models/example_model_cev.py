# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path
_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause
_install_double_click_pause()
import time
import numpy as np

import add_fp_to_path
import matplotlib.pyplot as plt

from financepy.models.cev import CEV
from financepy.utils.global_types import OptionTypes
from financepy.utils.global_vars import G_DAYS_IN_YEAR
from financepy.utils.date import Date



PLOT = False

########################################################################################


def test_analytical_models():

    value_dt = Date(1, 1, 2015)
    expiry_dt = Date(1, 1, 2016)

    interest_rate = 0.05
    dividend_yield = 0.02

    stock_price = 100.0
    spot_vol = 0.20

    tau = (expiry_dt - value_dt) / G_DAYS_IN_YEAR
    opt_type = OptionTypes.EUROPEAN_CALL.value

    num_steps = 100
    num_paths = 20000
    seed = 2838

    print(
        "TIME",
        "BETA",
        "SIGMA",
        "K",
        "FORMULA",
        "MC",
        "MCERR",
    )

    for beta in [0.25, 0.50, 0.75, 1.00]:

        # Choose sigma so that the initial percentage local volatility
        # is approximately spot_vol at S0.
        #
        # sigma_loc(S0) = sigma * S0^(beta - 1)
        #
        # therefore
        #
        # sigma = spot_vol * S0^(1 - beta)
        #
        sigma = spot_vol * stock_price ** (1.0 - beta)

        cev_model = CEV(
            sigma=sigma,
            beta=beta,
        )

        for strike_price in np.linspace(80.0, 120.0, 5):

            start = time.time()

            value_formula = cev_model.value(
                stock_price,
                tau,
                strike_price,
                opt_type,
                interest_rate,
                dividend_yield,
            )

            value_mc = cev_model.value_mc(
                stock_price,
                tau,
                strike_price,
                opt_type,
                interest_rate,
                dividend_yield,
                num_paths,
                num_steps,
                seed,
            )

            err = value_mc - value_formula

            end = time.time()
            elapsed = end - start

            print(
                f"{elapsed:6.3f}",
                f"{beta:7.5f}",
                f"{sigma:9.5f}",
                f"{strike_price:7.2f}",
                f"{value_formula:12.9f}",
                f"{value_mc:12.9f}",
                f"{err:12.9f}",
            )


########################################################################################


def test_monte_carlo():

    value_dt = Date(1, 1, 2015)
    expiry_dt = Date(1, 1, 2016)

    interest_rate = 0.05
    dividend_yield = 0.02

    stock_price = 100.0
    spot_vol = 0.20

    beta = 0.50

    sigma = spot_vol * stock_price ** (1.0 - beta)

    seed = 238

    tau = (expiry_dt - value_dt) / G_DAYS_IN_YEAR
    opt_type = OptionTypes.EUROPEAN_CALL.value

    cev_model = CEV(
        sigma=sigma,
        beta=beta,
    )

    print(
        "TIME",
        "BETA",
        "SIGMA",
        "K",
        "NSTEPS",
        "NPATHS",
        "FORMULA",
        "MCERR",
    )

    for strike_price in np.linspace(80.0, 120.0, 5):

        value_formula = cev_model.value(
            stock_price,
            tau,
            strike_price,
            opt_type,
            interest_rate,
            dividend_yield,
        )

        for num_steps in [25, 50, 100]:

            for num_paths in [10000, 20000]:

                start = time.time()

                value_mc = cev_model.value_mc(
                    stock_price,
                    tau,
                    strike_price,
                    opt_type,
                    interest_rate,
                    dividend_yield,
                    num_paths,
                    num_steps,
                    seed,
                )

                err_mc = value_mc - value_formula

                end = time.time()
                elapsed = end - start

                print(
                    elapsed,
                    beta,
                    sigma,
                    strike_price,
                    num_steps,
                    num_paths,
                    value_formula,
                    err_mc,
                )


########################################################################################


def test_cev_black_scholes_limit():
    """
    For beta = 1 the CEV model reduces to Black-Scholes.

    The implied-volatility curve should therefore be flat and equal
    to sigma.
    """

    model = CEV(
        sigma=0.20,
        beta=1.0,
    )

    t_exp = 1.0
    stock_price = 100.0
    interest_rate = 0.05
    dividend_yield = 0.02

    strikes = np.linspace(
        60.0,
        140.0,
        17,
    )

    vols = model.implied_volatility_curve(
        stock_price,
        t_exp,
        strikes,
        interest_rate,
        dividend_yield,
    )

    assert np.all(np.isfinite(vols))

    assert np.allclose(
        vols,
        0.20,
        rtol=1.0e-6,
        atol=1.0e-6,
    )

    print(
        "K",
        "IV",
    )

    for strike, vol in zip(strikes, vols):

        print(
            strike,
            vol,
        )


########################################################################################


def test_cev_volatility_skew():
    """
    Demonstrate the negative implied-volatility skew generated by
    beta < 1.
    """

    stock_price = 100.0
    t_exp = 1.0

    interest_rate = 0.05
    dividend_yield = 0.02

    spot_vol = 0.20
    beta = 0.50

    sigma = spot_vol * stock_price ** (1.0 - beta)

    model = CEV(
        sigma=sigma,
        beta=beta,
    )

    strikes = np.linspace(
        60.0,
        140.0,
        41,
    )

    vols = model.implied_volatility_curve(
        stock_price,
        t_exp,
        strikes,
        interest_rate,
        dividend_yield,
    )

    valid = np.isfinite(vols)

    assert np.all(valid)
    assert np.all(vols[valid] > 0.0)

    print(
        "K",
        "IV",
        "LOCAL_VOL_AT_K",
    )

    for strike, vol in zip(
        strikes,
        vols,
    ):

        local_vol = model.local_volatility(strike)

        print(
            strike,
            vol,
            local_vol,
        )

    if PLOT == True:

        plt.figure()

        plt.plot(
            strikes[valid],
            100.0 * vols[valid],
            marker="o",
            markersize=3,
        )

        plt.axvline(
            stock_price,
            linestyle="--",
        )

        plt.xlabel("Strike")
        plt.ylabel("Black-Scholes Implied Volatility (%)")

        plt.title("CEV Implied Volatility Skew")

        plt.grid(True)
        plt.show()


########################################################################################


def test_cev_beta_skews():
    """
    Compare the implied-volatility curves generated by different
    values of beta while keeping the initial local volatility fixed.
    """

    stock_price = 100.0
    t_exp = 1.0

    interest_rate = 0.05
    dividend_yield = 0.02

    spot_vol = 0.20

    strikes = np.linspace(
        60.0,
        140.0,
        41,
    )

    betas = [
        0.25,
        0.50,
        0.75,
        1.00,
    ]

    print(
        "BETA",
        "SIGMA",
        "K60",
        "K80",
        "K100",
        "K120",
        "K140",
    )

    if PLOT == True:
        plt.figure()

    for beta in betas:

        sigma = spot_vol * stock_price ** (1.0 - beta)

        model = CEV(
            sigma=sigma,
            beta=beta,
        )

        vols = model.implied_volatility_curve(
            stock_price,
            t_exp,
            strikes,
            interest_rate,
            dividend_yield,
        )

        assert np.all(np.isfinite(vols))

        assert np.all(vols > 0.0)

        selected_strikes = np.array(
            [
                60.0,
                80.0,
                100.0,
                120.0,
                140.0,
            ]
        )

        selected_vols = model.implied_volatility_curve(
            stock_price,
            t_exp,
            selected_strikes,
            interest_rate,
            dividend_yield,
        )

        print(
            beta,
            sigma,
            selected_vols[0],
            selected_vols[1],
            selected_vols[2],
            selected_vols[3],
            selected_vols[4],
        )

        if PLOT == True:

            plt.plot(
                strikes,
                100.0 * vols,
                marker="o",
                markersize=2,
                label=f"beta={beta:.2f}",
            )

    if PLOT == True:

        plt.axvline(
            stock_price,
            linestyle="--",
        )

        plt.xlabel("Strike")

        plt.ylabel("Black-Scholes Implied Volatility (%)")

        plt.title("CEV Implied Volatility Skew")

        plt.grid(True)
        plt.legend()
        plt.show()


########################################################################################


def test_cev_volatility_surface():

    stock_price = 100.0
    interest_rate = 0.05
    dividend_yield = 0.02

    spot_vol = 0.20
    beta = 0.50

    # Normalise sigma so that initial local volatility at S0 is 20%
    sigma = spot_vol * stock_price ** (1.0 - beta)

    model = CEV(
        sigma=sigma,
        beta=beta,
    )

    strikes = np.linspace(
        60.0,
        140.0,
        41,
    )

    expiries = np.array(
        [
            0.10,
            0.25,
            0.50,
            1.00,
            2.00,
            3.00,
            5.00,
        ]
    )

    vols = model.implied_volatility_surface(
        stock_price,
        expiries,
        strikes,
        interest_rate,
        dividend_yield,
    )

    assert vols.shape == (
        len(expiries),
        len(strikes),
    )

    assert np.all(np.isfinite(vols))

    assert np.all(vols > 0.0)

    print(
        "T",
        "K60",
        "K80",
        "K100",
        "K120",
        "K140",
    )

    selected_strikes = np.array(
        [
            60.0,
            80.0,
            100.0,
            120.0,
            140.0,
        ]
    )

    selected_vols = model.implied_volatility_surface(
        stock_price,
        expiries,
        selected_strikes,
        interest_rate,
        dividend_yield,
    )

    for i, t_exp in enumerate(expiries):

        print(
            t_exp,
            selected_vols[i, 0],
            selected_vols[i, 1],
            selected_vols[i, 2],
            selected_vols[i, 3],
            selected_vols[i, 4],
        )

    if PLOT == True:

        strike_grid, expiry_grid = np.meshgrid(
            strikes,
            expiries,
        )

        fig = plt.figure(figsize=(9, 6))

        ax = fig.add_subplot(
            111,
            projection="3d",
        )

        ax.plot_surface(
            strike_grid,
            expiry_grid,
            100.0 * vols,
        )

        ax.set_xlabel("Strike")
        ax.set_ylabel("Time to Expiry")
        ax.set_zlabel("Implied Volatility (%)")

        ax.set_title("CEV Implied Volatility Surface")

        plt.show()


########################################################################################


def test_cev_volatility_surfaces():

    stock_price = 100.0
    interest_rate = 0.05
    dividend_yield = 0.02

    spot_vol = 0.20

    strikes = np.linspace(
        60.0,
        140.0,
        41,
    )

    expiries = np.array(
        [
            0.25,
            0.50,
            1.00,
            2.00,
            3.00,
            5.00,
        ]
    )

    betas = [
        0.33333,
        0.66666,
        0.99999,
    ]

    strike_grid, expiry_grid = np.meshgrid(
        strikes,
        expiries,
    )

    if PLOT == True:

        fig = plt.figure(figsize=(10, 7))

        ax = fig.add_subplot(
            111,
            projection="3d",
        )

    for beta in betas:

        # Keep initial local volatility fixed at 20%
        sigma = spot_vol * stock_price ** (1.0 - beta)

        model = CEV(
            sigma=sigma,
            beta=beta,
        )

        vols = model.implied_volatility_surface(
            stock_price,
            expiries,
            strikes,
            interest_rate,
            dividend_yield,
        )

        assert vols.shape == (
            len(expiries),
            len(strikes),
        )

        assert np.all(np.isfinite(vols))

        assert np.all(vols > 0.0)

        print(
            "BETA",
            "T",
            "K60",
            "K80",
            "K100",
            "K120",
            "K140",
        )

        selected_strikes = np.array(
            [
                60.0,
                80.0,
                100.0,
                120.0,
                140.0,
            ]
        )

        selected_vols = model.implied_volatility_surface(
            stock_price,
            expiries,
            selected_strikes,
            interest_rate,
            dividend_yield,
        )

        for i, t_exp in enumerate(expiries):

            print(
                beta,
                t_exp,
                selected_vols[i, 0],
                selected_vols[i, 1],
                selected_vols[i, 2],
                selected_vols[i, 3],
                selected_vols[i, 4],
            )

        if PLOT == True:

            ax.plot_surface(
                strike_grid,
                expiry_grid,
                100.0 * vols,
                alpha=0.45,
                label=rf"$\beta={beta:.2f}$",
            )

    if PLOT == True:

        ax.set_xlabel("Strike")
        ax.set_ylabel("Time to Expiry")
        ax.set_zlabel("Implied Volatility (%)")

        ax.set_title("CEV Implied Volatility Surfaces")

        ax.view_init(
            elev=25,
            azim=-130,
        )

        plt.show()


########################################################################################


test_analytical_models()
test_monte_carlo()
test_cev_black_scholes_limit()
test_cev_volatility_skew()
test_cev_beta_skews()
test_cev_volatility_surface()
test_cev_volatility_surfaces()

