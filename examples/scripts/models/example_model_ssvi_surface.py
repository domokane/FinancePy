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

from financepy.models.ssvi_surface import SSVIPowerLawPhi
from financepy.models.ssvi_surface import SSVISurface



########################################################################################
########################################################################################
# TEST SSVI SURFACE
########################################################################################
########################################################################################


def ssvi_total_variance(
    k,
    theta,
    rho,
    eta,
    gamma,
):

    phi = eta * theta ** (-gamma)

    x = phi * k

    return 0.5 * theta * (1.0 + rho * x + np.sqrt((x + rho) ** 2 + 1.0 - rho * rho))


########################################################################################
# TEST POWER-LAW PHI
########################################################################################


def test_ssvi_power_law_phi():

    eta = 1.20
    gamma = 0.40

    phi = SSVIPowerLawPhi(
        eta,
        gamma,
    )

    thetas = np.array(
        [
            0.01,
            0.02,
            0.05,
            0.10,
            0.20,
        ]
    )

    print(
        "THETA",
        "EXPECTED",
        "CALCULATED",
        "ERROR",
    )

    for theta in thetas:

        expected = eta * theta ** (-gamma)

        calculated = phi(theta)

        error = calculated - expected

        print(
            theta,
            expected,
            calculated,
            error,
        )

        assert abs(error) < 1.0e-12


########################################################################################
# TEST SSVI TOTAL VARIANCE FORMULA
########################################################################################


def test_ssvi_formula():

    expiries = np.array(
        [
            1.00,
        ]
    )

    theta = np.array(
        [
            0.04,
        ]
    )

    rho = -0.70
    eta = 1.20
    gamma = 0.40

    phi = SSVIPowerLawPhi(
        eta,
        gamma,
    )

    surface = SSVISurface(
        expiries,
        theta,
        rho,
        phi,
    )

    ks = np.array(
        [
            -0.40,
            -0.30,
            -0.20,
            -0.10,
            0.00,
            0.10,
            0.20,
            0.30,
            0.40,
        ]
    )

    print(
        "K",
        "EXPECTED W",
        "SSVI W",
        "ERROR",
    )

    for k in ks:

        expected = ssvi_total_variance(
            k,
            theta[0],
            rho,
            eta,
            gamma,
        )

        calculated = surface.total_variance_from_log_moneyness(
            k,
            expiries[0],
        )

        error = calculated - expected

        print(
            k,
            expected,
            calculated,
            error,
        )

        assert abs(error) < 1.0e-12


########################################################################################
# TEST ATM PROPERTY
########################################################################################


def test_ssvi_atm():

    expiries = np.array(
        [
            0.50,
            1.00,
            2.00,
            5.00,
        ]
    )

    theta = np.array(
        [
            0.020,
            0.040,
            0.080,
            0.180,
        ]
    )

    phi = SSVIPowerLawPhi(
        eta=1.00,
        gamma=0.30,
    )

    surface = SSVISurface(
        expiries,
        theta,
        rho=-0.60,
        phi_function=phi,
    )

    print(
        "T",
        "THETA",
        "W(0,T)",
        "ERROR",
    )

    for i, t_exp in enumerate(expiries):

        calculated = surface.total_variance_from_log_moneyness(
            0.0,
            t_exp,
        )

        expected = theta[i]

        error = calculated - expected

        print(
            t_exp,
            expected,
            calculated,
            error,
        )

        assert abs(error) < 1.0e-12


########################################################################################
# TEST FORWARD SCALE INVARIANCE
########################################################################################


def test_ssvi_forward_scale_invariance():

    expiries = np.array(
        [
            1.00,
        ]
    )

    theta = np.array(
        [
            0.04,
        ]
    )

    phi = SSVIPowerLawPhi(
        eta=1.20,
        gamma=0.40,
    )

    surface = SSVISurface(
        expiries,
        theta,
        rho=-0.70,
        phi_function=phi,
    )

    t_exp = 1.00

    vol_1 = surface.implied_volatility(
        100.0,
        90.0,
        t_exp,
    )

    vol_2 = surface.implied_volatility(
        200.0,
        180.0,
        t_exp,
    )

    vol_3 = surface.implied_volatility(
        50.0,
        45.0,
        t_exp,
    )

    print(
        "VOL 1",
        "VOL 2",
        "VOL 3",
    )

    print(
        vol_1,
        vol_2,
        vol_3,
    )

    assert abs(vol_1 - vol_2) < 1.0e-12

    assert abs(vol_1 - vol_3) < 1.0e-12


########################################################################################
# TEST THETA INTERPOLATION
########################################################################################


def test_ssvi_theta_interpolation():

    expiries = np.array(
        [
            1.00,
            2.00,
        ]
    )

    theta = np.array(
        [
            0.04,
            0.10,
        ]
    )

    phi = SSVIPowerLawPhi(
        eta=1.10,
        gamma=0.30,
    )

    surface = SSVISurface(
        expiries,
        theta,
        rho=-0.50,
        phi_function=phi,
    )

    t_exp = 1.50

    expected = 0.07

    calculated = surface.theta(t_exp)

    error = calculated - expected

    print(
        "T",
        "EXPECTED THETA",
        "SSVI THETA",
        "ERROR",
    )

    print(
        t_exp,
        expected,
        calculated,
        error,
    )

    assert abs(error) < 1.0e-12


########################################################################################
# TEST INTERPOLATED SSVI SLICE
########################################################################################


def test_ssvi_surface_interpolation():

    expiries = np.array(
        [
            1.00,
            2.00,
        ]
    )

    theta = np.array(
        [
            0.040,
            0.090,
        ]
    )

    rho = -0.60
    eta = 1.10
    gamma = 0.35

    phi = SSVIPowerLawPhi(
        eta,
        gamma,
    )

    surface = SSVISurface(
        expiries,
        theta,
        rho,
        phi,
    )

    forward = 100.0
    strike = 90.0
    t_exp = 1.50

    k = np.log(strike / forward)

    expected_theta = theta[0] + 0.5 * (theta[1] - theta[0])

    expected_w = ssvi_total_variance(
        k,
        expected_theta,
        rho,
        eta,
        gamma,
    )

    expected_vol = np.sqrt(expected_w / t_exp)

    calculated_w = surface.total_variance(
        forward,
        strike,
        t_exp,
    )

    calculated_vol = surface.implied_volatility(
        forward,
        strike,
        t_exp,
    )

    print(
        "EXPECTED W",
        "SSVI W",
        "EXPECTED VOL",
        "SSVI VOL",
    )

    print(
        expected_w,
        calculated_w,
        expected_vol,
        calculated_vol,
    )

    assert abs(calculated_w - expected_w) < 1.0e-12

    assert abs(calculated_vol - expected_vol) < 1.0e-12


########################################################################################
# TEST SURFACE GRID
########################################################################################


def test_ssvi_surface_grid():

    expiries = np.array(
        [
            0.50,
            1.00,
            2.00,
        ]
    )

    theta = np.array(
        [
            0.020,
            0.040,
            0.090,
        ]
    )

    phi = SSVIPowerLawPhi(
        eta=1.10,
        gamma=0.35,
    )

    surface = SSVISurface(
        expiries,
        theta,
        rho=-0.60,
        phi_function=phi,
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

    print(
        "NUM EXPIRIES",
        "NUM STRIKES",
        "MIN VOL",
        "MAX VOL",
    )

    print(
        vols.shape[0],
        vols.shape[1],
        np.min(vols),
        np.max(vols),
    )

    assert vols.shape == (
        len(query_expiries),
        len(strikes),
    )

    assert np.all(np.isfinite(vols))

    assert np.all(vols > 0.0)


########################################################################################
# TEST PARAMETER ACCESSORS
########################################################################################


def test_ssvi_parameters():

    expiries = np.array(
        [
            0.50,
            1.00,
            2.00,
        ]
    )

    theta = np.array(
        [
            0.020,
            0.040,
            0.080,
        ]
    )

    rho = -0.65
    eta = 1.15
    gamma = 0.35

    phi = SSVIPowerLawPhi(
        eta,
        gamma,
    )

    surface = SSVISurface(
        expiries,
        theta,
        rho,
        phi,
    )

    parameters = surface.parameters()

    returned_theta = surface.atm_total_variances()

    returned_expiries = surface.expiries()

    print(
        "RHO",
        "ETA",
        "GAMMA",
    )

    print(
        parameters[0],
        parameters[1],
        parameters[2],
    )

    assert abs(parameters[0] - rho) < 1.0e-12

    assert abs(parameters[1] - eta) < 1.0e-12

    assert abs(parameters[2] - gamma) < 1.0e-12

    assert np.max(np.abs(returned_theta - theta)) < 1.0e-12

    assert np.max(np.abs(returned_expiries - expiries)) < 1.0e-12


########################################################################################
# TEST JOINT CALIBRATION
########################################################################################


def test_ssvi_surface_calibration():

    ####################################################################################
    # KNOWN SYNTHETIC SSVI SURFACE
    ####################################################################################

    expiries = np.array(
        [
            0.50,
            1.00,
            2.00,
            5.00,
        ]
    )

    forwards = np.array(
        [
            100.0,
            100.0,
            100.0,
            100.0,
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

    expected_theta = np.array(
        [
            0.020,
            0.040,
            0.080,
            0.180,
        ]
    )

    expected_rho = -0.60
    expected_eta = 1.00
    expected_gamma = 0.30

    phi = SSVIPowerLawPhi(
        expected_eta,
        expected_gamma,
    )

    expected_surface = SSVISurface(
        expiries,
        expected_theta,
        expected_rho,
        phi,
    )

    ####################################################################################
    # GENERATE EXACT SYNTHETIC MARKET VOLATILITIES
    ####################################################################################

    market_vols = expected_surface.implied_volatility_surface(
        forwards,
        strikes,
        expiries,
    )

    ####################################################################################
    # CALIBRATE A NEW SSVI SURFACE
    ####################################################################################

    calibrated_surface = SSVISurface()

    cost = calibrated_surface.calibrate(
        forwards,
        strikes,
        expiries,
        market_vols,
    )

    calibrated_theta = calibrated_surface.atm_total_variances()

    (
        calibrated_rho,
        calibrated_eta,
        calibrated_gamma,
    ) = calibrated_surface.parameters()

    ####################################################################################
    # CALCULATE FITTED SURFACE
    ####################################################################################

    fitted_vols = calibrated_surface.implied_volatility_surface(
        forwards,
        strikes,
        expiries,
    )

    max_vol_error = np.max(np.abs(fitted_vols - market_vols))

    max_theta_error = np.max(np.abs(calibrated_theta - expected_theta))

    ####################################################################################
    # OUTPUT
    ####################################################################################

    print(
        "PARAMETER",
        "EXPECTED",
        "CALIBRATED",
        "ERROR",
    )

    print(
        "RHO",
        expected_rho,
        calibrated_rho,
        calibrated_rho - expected_rho,
    )

    print(
        "ETA",
        expected_eta,
        calibrated_eta,
        calibrated_eta - expected_eta,
    )

    print(
        "GAMMA",
        expected_gamma,
        calibrated_gamma,
        calibrated_gamma - expected_gamma,
    )

    print(
        "T",
        "EXPECTED THETA",
        "CALIBRATED THETA",
        "ERROR",
    )

    for i, t_exp in enumerate(expiries):

        print(
            t_exp,
            expected_theta[i],
            calibrated_theta[i],
            calibrated_theta[i] - expected_theta[i],
        )

    print(
        "COST",
        "MAX THETA ERROR",
        "MAX VOL ERROR",
    )

    print(
        cost,
        max_theta_error,
        max_vol_error,
    )

    ####################################################################################
    # ASSERTIONS
    ####################################################################################

    assert np.isfinite(cost)

    assert max_vol_error < 1.0e-7

    assert max_theta_error < 1.0e-6

    assert abs(calibrated_rho - expected_rho) < 1.0e-5

    assert abs(calibrated_eta - expected_eta) < 1.0e-5

    assert abs(calibrated_gamma - expected_gamma) < 1.0e-5


########################################################################################
# TEST CALIBRATED THETA MONOTONICITY
########################################################################################


def test_ssvi_calibrated_theta_monotonicity():

    expiries = np.array(
        [
            0.25,
            0.50,
            1.00,
            2.00,
            5.00,
        ]
    )

    forwards = np.array(
        [
            100.0,
            100.0,
            100.0,
            100.0,
            100.0,
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

    market_vols = np.array(
        [
            [0.300, 0.270, 0.240, 0.210, 0.195, 0.187, 0.185],
            [0.285, 0.260, 0.232, 0.205, 0.192, 0.185, 0.184],
            [0.270, 0.250, 0.225, 0.200, 0.190, 0.185, 0.185],
            [0.250, 0.235, 0.215, 0.198, 0.191, 0.189, 0.190],
            [0.230, 0.220, 0.208, 0.197, 0.193, 0.191, 0.192],
        ]
    )

    surface = SSVISurface()

    surface.calibrate(
        forwards,
        strikes,
        expiries,
        market_vols,
    )

    theta = surface.atm_total_variances()

    increments = np.diff(theta)

    print(
        "T",
        "THETA",
    )

    for t_exp, theta_value in zip(
        expiries,
        theta,
    ):

        print(
            t_exp,
            theta_value,
        )

    assert np.all(theta > 0.0)

    assert np.all(increments > 0.0)


########################################################################################
# RUN TESTS
########################################################################################


test_ssvi_power_law_phi()
test_ssvi_formula()
test_ssvi_atm()
test_ssvi_forward_scale_invariance()
test_ssvi_theta_interpolation()
test_ssvi_surface_interpolation()
test_ssvi_surface_grid()
test_ssvi_parameters()
test_ssvi_surface_calibration()
test_ssvi_calibrated_theta_monotonicity()

