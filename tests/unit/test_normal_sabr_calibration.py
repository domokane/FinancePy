import numpy as np
import pytest
from financepy.models.sabr_normal import NormalSABR, normal_sabr_vol
from financepy.utils.error import FinError


def test_normal_sabr_flat_normal_vol_for_zero_vol_of_vol():
    for strike in (-0.03, 0.0, 0.01, 0.03):
        assert normal_sabr_vol(-0.005, strike, 3.0, 0.012, -0.3, 0.0) == pytest.approx(0.012)


def test_normal_sabr_atm_limit_and_negative_rates():
    alpha, rho, nu, t = 0.015, -0.4, 0.6, 2.0
    expected = alpha * (1 + (2 - 3 * rho**2) * nu**2 * t / 24)
    assert normal_sabr_vol(-0.012, -0.012, t, alpha, rho, nu) == pytest.approx(expected)
    assert normal_sabr_vol(-0.012, -0.012 + 1e-11, t, alpha, rho, nu) == pytest.approx(expected, rel=1e-7)


@pytest.mark.parametrize("rho", [-0.8, -0.2, 0.0, 0.7])
def test_normal_sabr_is_positive_and_continuous(rho):
    strikes = np.linspace(-0.05, 0.05, 41)
    vols = np.array([normal_sabr_vol(-0.01, k, 1.5, 0.012, rho, 0.8) for k in strikes])
    assert np.all(np.isfinite(vols)) and np.all(vols > 0)


def test_normal_sabr_calibration_recovers_synthetic_smile():
    true = NormalSABR(0.012, -0.35, 0.70)
    strikes = np.linspace(-0.04, 0.04, 13)
    market = [true.normal_vol(-0.008, float(k), 2.0) for k in strikes]
    calibrated = NormalSABR.calibrate(
        -0.008, strikes, market, 2.0, initial_guess=(0.011, -0.2, 0.4)
    )
    np.testing.assert_allclose(
        [calibrated.alpha, calibrated.rho, calibrated.nu],
        [true.alpha, true.rho, true.nu], rtol=2e-5, atol=1e-7
    )


@pytest.mark.parametrize("alpha,rho,nu", [(0,0,0.2), (0.01,1.0,0.2), (0.01,0,-0.2)])
def test_normal_sabr_rejects_bad_params(alpha,rho,nu):
    with pytest.raises(FinError):
        NormalSABR(alpha,rho,nu)


def test_calibration_rejects_too_few_quotes():
    with pytest.raises(FinError):
        NormalSABR.calibrate(0.01, [0.0, 0.02], [0.01, 0.01], 1.0)

