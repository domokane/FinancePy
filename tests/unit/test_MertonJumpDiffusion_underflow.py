import pytest

from financepy.models.black_scholes_analytic import value as black_scholes_value
from financepy.models.merton_jump_diffusion import MertonJumpDiffusion
from financepy.utils.global_types import OptionTypes


def test_zero_jump_size_preserves_price_at_large_poisson_exposure():
    model = MertonJumpDiffusion(
        sigma=0.2,
        jump_intensity=800.0,
        jump_mean=0.0,
        jump_volatility=0.0,
        max_jumps=4000,
    )

    price = model.value(120.0, 1.0, 100.0, 0.03, 0.01, OptionTypes.EUROPEAN_CALL)
    expected = black_scholes_value(
        120.0,
        1.0,
        100.0,
        0.03,
        0.01,
        0.2,
        OptionTypes.EUROPEAN_CALL.value,
    )

    assert price == pytest.approx(expected, rel=1.0e-9)
