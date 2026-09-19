
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

from financepy.utils.global_types import OptionTypes
from financepy.models.black import Black
from financepy.models.lognormal_mixture_model import LognormalMixtureModel
from financepy.models.lognormal_mixture_surface import LognormalMixtureSurface




def test_forward_constraint():
    model = LognormalMixtureModel(F=100.0, T=1.0, r=0.03)

    p = 0.4
    displacement = -0.05

    F1, F2 = model._component_forwards(p, displacement)

    mixture_forward = p * F1 + (1.0 - p) * F2

    assert np.isclose(mixture_forward, model.F, atol=1e-12)


def test_mixture_reduces_to_black():
    """
    If both mixture components have the same forward and volatility,
    the mixture price must equal the FinancePy Black price.
    """

    F = 100.0
    T = 1.0
    r = 0.03
    sigma = 0.20

    strikes = np.array([80.0, 90.0, 100.0, 110.0, 120.0])

    model = LognormalMixtureModel(F=F, T=T, r=r)

    mixture_prices = model.price(K=strikes, p=0.35, displacement=0.0, sigma1=sigma, sigma2=sigma)

    black = Black(sigma)

    discount_factor = np.exp(-r * T)

    black_prices = np.array([black.value(F, K, T, discount_factor, OptionTypes.EUROPEAN_CALL) for K in strikes])

    assert np.allclose(mixture_prices, black_prices, atol=1e-10)


def test_mixture_price_decreases_with_strike():
    model = LognormalMixtureModel(F=100.0, T=1.0, r=0.03)

    strikes = np.linspace(70.0, 130.0, 50)

    prices = model.price(K=strikes, p=0.4, displacement=-0.05, sigma1=0.15, sigma2=0.30)

    assert np.all(np.diff(prices) <= 1e-12)


def test_mixture_price_convex_in_strike():
    model = LognormalMixtureModel(F=100.0, T=1.0, r=0.03)

    strikes = np.linspace(70.0, 130.0, 100)

    prices = model.price(K=strikes, p=0.4, displacement=-0.05, sigma1=0.15, sigma2=0.30)

    second_diff = prices[:-2] - 2.0 * prices[1:-1] + prices[2:]

    assert np.all(second_diff >= -1e-10)


def test_density_positive():
    model = LognormalMixtureModel(F=100.0, T=1.0)

    model.p = 0.4
    model.displacement = -0.05
    model.sigma1 = 0.15
    model.sigma2 = 0.30

    model.F1, model.F2 = model._component_forwards(model.p, model.displacement)

    s = np.linspace(1e-4, 400.0, 10000)

    density = model.density(s)

    assert np.all(density >= 0.0)


def test_density_integrates_to_one():
    model = LognormalMixtureModel(F=100.0, T=1.0)

    model.p = 0.4
    model.displacement = -0.05
    model.sigma1 = 0.15
    model.sigma2 = 0.30

    model.F1, model.F2 = model._component_forwards(model.p, model.displacement)

    s = np.linspace(1e-4, 500.0, 100000)

    density = model.density(s)

    integral = np.trapezoid(density, s)

    assert np.isclose(integral, 1.0, atol=1e-4)


def test_density_mean_equals_forward():
    model = LognormalMixtureModel(F=100.0, T=1.0)

    model.p = 0.4
    model.displacement = -0.05
    model.sigma1 = 0.15
    model.sigma2 = 0.30

    model.F1, model.F2 = model._component_forwards(model.p, model.displacement)

    s = np.linspace(1e-4, 500.0, 100000)

    density = model.density(s)

    expected_s = np.trapezoid(s * density, s)

    assert np.isclose(expected_s, model.F, atol=1e-3)


def test_calibration_recovers_synthetic_smile():
    """
    Generate a synthetic smile from known mixture parameters
    and check that calibration reproduces the smile.
    """

    F = 100.0
    T = 1.0
    r = 0.03

    true_model = LognormalMixtureModel(F=F, T=T, r=r)

    true_model.p = 0.35
    true_model.displacement = -0.06
    true_model.sigma1 = 0.14
    true_model.sigma2 = 0.28

    true_model.F1, true_model.F2 = true_model._component_forwards(true_model.p, true_model.displacement)

    strikes = np.linspace(75.0, 130.0, 20)

    market_vols = true_model.implied_vol(strikes)

    calibrated_model = LognormalMixtureModel(F=F, T=T, r=r)

    calibrated_model.calibrate(strikes, market_vols)

    fitted_vols = calibrated_model.implied_vol(strikes)

    assert np.max(np.abs(fitted_vols - market_vols)) < 1e-4


#########################################################################################


def test_surface_calibration():

    # Replace these with the actual surface used earlier
    maturities = np.array([0.25, 0.50, 1.00, 2.00])

    forwards = np.array([100.5, 101.0, 102.0, 104.0])

    rates = np.array([0.03, 0.03, 0.032, 0.035])

    strikes = np.array([80.0, 90.0, 100.0, 110.0, 120.0])

    # Actual volatility surface
    vol_surface = np.array(
        [
            [0.3000, 0.2600, 0.2300, 0.2200, 0.2300],
            [0.2800, 0.2400, 0.2100, 0.2000, 0.2100],
            [0.2500, 0.2200, 0.2000, 0.1900, 0.2000],
            [0.2300, 0.2100, 0.1900, 0.1900, 0.2000],
        ]
    )

    # --------------------------------------------------------
    # Calibrate to the volatility surface
    # --------------------------------------------------------

    surface = LognormalMixtureSurface(
        maturities=maturities, forwards=forwards, strikes=strikes, vol_surface=vol_surface, rates=rates
    )

    surface.calibrate()

    # --------------------------------------------------------
    # Reconstruct fitted surface
    # --------------------------------------------------------

    fitted_surface = surface.fitted_vol_surface()

    errors = fitted_surface - vol_surface

    max_error = np.max(np.abs(errors))

    rmse = np.sqrt(np.mean(errors**2))

    # --------------------------------------------------------
    # Regression checks
    # --------------------------------------------------------

    assert np.all(np.isfinite(fitted_surface))

    assert np.all(fitted_surface > 0.0)

    # Reasonable calibration quality
    assert max_error < 0.01
    assert rmse < 0.005

    # Forward condition at each maturity
    for i, model in enumerate(surface.models):

        assert np.isclose(model.forward_check, forwards[i], atol=1e-10)

    # --------------------------------------------------------
    # Diagnostics
    # --------------------------------------------------------

    if 1 == 0:
        print("Market surface:")
        print(vol_surface)

        print("Fitted surface:")
        print(fitted_surface)

        print("Errors:")
        print(errors)

        print("Maximum calibration error:", max_error)
        print("Surface RMSE:", rmse)

        for i, model in enumerate(surface.models):

            print(maturities[i], model.p, model.F1, model.F2, model.sigma1, model.sigma2)


#########################################################################################

test_forward_constraint()
test_mixture_reduces_to_black()
test_mixture_price_decreases_with_strike()
test_mixture_price_convex_in_strike()
test_density_positive()
test_density_integrates_to_one()
test_density_mean_equals_forward()
test_calibration_recovers_synthetic_smile()
test_surface_calibration()

