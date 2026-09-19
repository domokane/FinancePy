########################################################################################
########################################################################################
# SVI REGRESSION TESTS
########################################################################################
########################################################################################


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

from financepy.models.svi import SVI
from financepy.models.svi_surface import SVISurface




########################################################################################
# RAW SVI FORMULA
########################################################################################


def raw_svi(
    k,
    a,
    b,
    rho,
    m,
    sigma,
):

    x = k - m

    return a + b * (rho * x + np.sqrt(x * x + sigma * sigma))


########################################################################################
# SVI FORMULA
########################################################################################


def test_svi_formula():

    a = 0.020
    b = 0.100
    rho = -0.40
    m = 0.050
    sigma = 0.200

    svi = SVI(
        a,
        b,
        rho,
        m,
        sigma,
    )

    forward = 100.0
    t_exp = 1.50

    log_moneyness = np.array(
        [
            -0.30,
            -0.20,
            -0.10,
            0.00,
            0.10,
            0.20,
            0.30,
        ]
    )

    strikes = forward * np.exp(log_moneyness)

    print(
        "K",
        "LOG-MONEYNESS",
        "EXPECTED W",
        "SVI W",
        "EXPECTED VOL",
        "SVI VOL",
    )

    for k, strike in zip(
        log_moneyness,
        strikes,
    ):

        expected_w = raw_svi(
            k,
            a,
            b,
            rho,
            m,
            sigma,
        )

        svi_w = svi.total_variance(
            forward,
            strike,
        )

        expected_vol = np.sqrt(expected_w / t_exp)

        svi_vol = svi.implied_volatility(
            forward,
            strike,
            t_exp,
        )

        print(
            strike,
            k,
            expected_w,
            svi_w,
            expected_vol,
            svi_vol,
        )

        assert abs(svi_w - expected_w) < 1.0e-12

        assert abs(svi_vol - expected_vol) < 1.0e-12


########################################################################################
# SVI LOG-MONEYNESS
########################################################################################


def test_svi_log_moneyness():

    a = 0.020
    b = 0.100
    rho = -0.40
    m = 0.050
    sigma = 0.200

    svi = SVI(
        a,
        b,
        rho,
        m,
        sigma,
    )

    ks = np.linspace(
        -0.50,
        0.50,
        21,
    )

    for k in ks:

        expected_w = raw_svi(
            k,
            a,
            b,
            rho,
            m,
            sigma,
        )

        svi_w = svi.total_variance_from_log_moneyness(k)

        assert abs(svi_w - expected_w) < 1.0e-12


########################################################################################
# FORWARD SCALE INVARIANCE
########################################################################################


def test_svi_forward_scale_invariance():

    svi = SVI(
        a=0.020,
        b=0.100,
        rho=-0.40,
        m=0.050,
        sigma=0.200,
    )

    t_exp = 2.0

    vol_1 = svi.implied_volatility(
        100.0,
        90.0,
        t_exp,
    )

    vol_2 = svi.implied_volatility(
        200.0,
        180.0,
        t_exp,
    )

    vol_3 = svi.implied_volatility(
        50.0,
        45.0,
        t_exp,
    )

    assert abs(vol_1 - vol_2) < 1.0e-12

    assert abs(vol_1 - vol_3) < 1.0e-12


########################################################################################
# SVI CALIBRATION
########################################################################################


def test_svi_calibration():

    forward = 100.0
    t_exp = 1.0

    expected_parameters = np.array(
        [
            0.020,
            0.100,
            -0.40,
            0.050,
            0.200,
        ]
    )

    (
        a,
        b,
        rho,
        m,
        sigma,
    ) = expected_parameters

    strikes = np.linspace(
        70.0,
        130.0,
        21,
    )

    market_vols = np.empty(
        len(strikes),
        dtype=float,
    )

    for i, strike in enumerate(strikes):

        k = np.log(strike / forward)

        w = raw_svi(
            k,
            a,
            b,
            rho,
            m,
            sigma,
        )

        market_vols[i] = np.sqrt(w / t_exp)

    svi = SVI()

    svi.calibrate(
        forward,
        strikes,
        market_vols,
        t_exp,
    )

    calibrated_parameters = svi.parameters()

    print(
        "PARAMETER",
        "EXPECTED",
        "CALIBRATED",
        "ERROR",
    )

    names = [
        "A",
        "B",
        "RHO",
        "M",
        "SIGMA",
    ]

    for name, expected, calibrated in zip(
        names,
        expected_parameters,
        calibrated_parameters,
    ):

        error = calibrated - expected

        print(
            name,
            expected,
            calibrated,
            error,
        )

    #
    # Parameter recovery can be less robust than price/smile recovery
    # because SVI parameters can be highly correlated. The stronger
    # regression test is therefore reproduction of the generated smile.
    #

    fitted_vols = svi.implied_volatility_curve(
        forward,
        strikes,
        t_exp,
    )

    max_error = np.max(np.abs(fitted_vols - market_vols))

    assert max_error < 1.0e-6


########################################################################################
# SVI SURFACE CALIBRATION
########################################################################################


def test_svi_surface_calibration():

    expiries = np.array(
        [
            0.50,
            1.00,
            2.00,
        ]
    )

    forwards = np.array(
        [
            100.0,
            102.0,
            105.0,
        ]
    )

    strikes = np.array(
        [
            70.0,
            80.0,
            90.0,
            100.0,
            110.0,
            120.0,
            130.0,
        ]
    )

    parameters = np.array(
        [
            [
                0.010,
                0.060,
                -0.30,
                0.020,
                0.150,
            ],
            [
                0.020,
                0.100,
                -0.40,
                0.050,
                0.200,
            ],
            [
                0.030,
                0.120,
                -0.50,
                0.070,
                0.250,
            ],
        ]
    )

    market_vols = np.empty(
        (
            len(expiries),
            len(strikes),
        ),
        dtype=float,
    )

    ####################################################################################
    # GENERATE SYNTHETIC MARKET SURFACE
    ####################################################################################

    for i, t_exp in enumerate(expiries):

        forward = forwards[i]

        (
            a,
            b,
            rho,
            m,
            sigma,
        ) = parameters[i]

        for j, strike in enumerate(strikes):

            k = np.log(strike / forward)

            w = raw_svi(
                k,
                a,
                b,
                rho,
                m,
                sigma,
            )

            market_vols[i, j] = np.sqrt(w / t_exp)

    ####################################################################################
    # CALIBRATE SURFACE
    ####################################################################################

    surface = SVISurface()

    calibration_errors = surface.calibrate(
        forwards,
        strikes,
        expiries,
        market_vols,
    )

    ####################################################################################
    # CHECK CALIBRATION NODES
    ####################################################################################

    print(
        "T",
        "K",
        "MARKET VOL",
        "SVI VOL",
        "ERROR",
    )

    max_error = 0.0

    for i, t_exp in enumerate(expiries):

        forward = forwards[i]

        for j, strike in enumerate(strikes):

            svi_vol = surface.implied_volatility(
                forward,
                strike,
                t_exp,
            )

            market_vol = market_vols[
                i,
                j,
            ]

            error = svi_vol - market_vol

            max_error = max(
                max_error,
                abs(error),
            )

            print(
                t_exp,
                strike,
                market_vol,
                svi_vol,
                error,
            )

    assert max_error < 1.0e-6

    assert np.all(np.isfinite(calibration_errors))


########################################################################################
# SVI SURFACE INTERPOLATION
########################################################################################


def test_svi_surface_interpolation():

    expiries = np.array(
        [
            1.00,
            2.00,
        ]
    )

    parameters = np.array(
        [
            [
                0.020,
                0.100,
                -0.40,
                0.050,
                0.200,
            ],
            [
                0.030,
                0.120,
                -0.50,
                0.070,
                0.230,
            ],
        ]
    )

    surface = SVISurface(
        expiries,
        parameters,
    )

    forward = 100.0
    strike = 90.0

    t_exp = 1.50

    k = np.log(strike / forward)

    ####################################################################################
    # LOWER SLICE
    ####################################################################################

    w1 = raw_svi(
        k,
        *parameters[0],
    )

    ####################################################################################
    # UPPER SLICE
    ####################################################################################

    w2 = raw_svi(
        k,
        *parameters[1],
    )

    ####################################################################################
    # EXPECTED LINEAR INTERPOLATION IN TOTAL VARIANCE
    ####################################################################################

    alpha = (t_exp - expiries[0]) / (expiries[1] - expiries[0])

    expected_w = (1.0 - alpha) * w1 + alpha * w2

    surface_w = surface.total_variance(
        forward,
        strike,
        t_exp,
    )

    expected_vol = np.sqrt(expected_w / t_exp)

    surface_vol = surface.implied_volatility(
        forward,
        strike,
        t_exp,
    )

    print(
        "EXPECTED W",
        "SURFACE W",
        "EXPECTED VOL",
        "SURFACE VOL",
    )

    print(
        expected_w,
        surface_w,
        expected_vol,
        surface_vol,
    )

    assert abs(surface_w - expected_w) < 1.0e-12

    assert abs(surface_vol - expected_vol) < 1.0e-12


########################################################################################
# SVI SURFACE GRID
########################################################################################


def test_svi_surface_grid():

    expiries = np.array(
        [
            0.50,
            1.00,
            2.00,
        ]
    )

    parameters = np.array(
        [
            [
                0.010,
                0.060,
                -0.30,
                0.020,
                0.150,
            ],
            [
                0.020,
                0.100,
                -0.40,
                0.050,
                0.200,
            ],
            [
                0.030,
                0.120,
                -0.50,
                0.070,
                0.250,
            ],
        ]
    )

    surface = SVISurface(
        expiries,
        parameters,
    )

    query_expiries = np.array(
        [
            0.50,
            0.75,
            1.00,
            1.50,
            2.00,
        ]
    )

    forwards = np.array(
        [
            100.0,
            101.0,
            102.0,
            103.0,
            105.0,
        ]
    )

    strikes = np.array(
        [
            80.0,
            90.0,
            100.0,
            110.0,
            120.0,
        ]
    )

    vols = surface.implied_volatility_surface(
        forwards,
        strikes,
        query_expiries,
    )

    assert vols.shape == (
        len(query_expiries),
        len(strikes),
    )

    assert np.all(np.isfinite(vols))

    assert np.all(vols > 0.0)


########################################################################################
# RUN TESTS
########################################################################################


test_svi_formula()
test_svi_log_moneyness()
test_svi_forward_scale_invariance()
test_svi_calibration()
test_svi_surface_calibration()
test_svi_surface_interpolation()
test_svi_surface_grid()


