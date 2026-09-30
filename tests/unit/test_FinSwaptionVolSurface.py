import numpy as np
import pytest

from financepy.utils.date import Date
from financepy.market.volatility.swaption_vol_surface import SwaptionVolSurface
from financepy.utils.global_types import VolFuncTypes


def _surface(vol_func_type=VolFuncTypes.SABR_BETA_HALF):
    value_dt = Date(1, 1, 2026)
    expiry_dts = [value_dt.add_years(y) for y in (1, 2, 5)]
    fwd_swap_rates = np.array([0.030, 0.032, 0.035])
    offsets = np.array([-0.010, -0.005, 0.0, 0.005, 0.010])

    # First dimension is the strike, then the expiry date
    strike_grid = offsets[:, None] + fwd_swap_rates[None, :]
    vol_grid = np.array(
        [
            [0.300, 0.290, 0.280],
            [0.270, 0.265, 0.260],
            [0.250, 0.245, 0.240],
            [0.260, 0.255, 0.250],
            [0.280, 0.275, 0.270],
        ]
    )

    surface = SwaptionVolSurface(
        value_dt,
        expiry_dts,
        fwd_swap_rates,
        strike_grid,
        vol_grid,
        vol_func_type,
    )
    return surface, expiry_dts, strike_grid, vol_grid


@pytest.mark.parametrize(
    "vol_func_type",
    [VolFuncTypes.SABR, VolFuncTypes.SABR_BETA_HALF, VolFuncTypes.SABR_BETA_ONE],
)
def test_check_calibration_accepts_numpy_grids(vol_func_type):
    surface, expiry_dts, strike_grid, vol_grid = _surface(vol_func_type)

    max_abs_diff = surface.check_calibration(False)

    expected = max(
        abs(surface.vol_from_strike_dt(strike_grid[j][i], expiry_dts[i]) - vol_grid[j][i])
        for i in range(len(expiry_dts))
        for j in range(len(strike_grid))
    )
    assert max_abs_diff == expected
    assert max_abs_diff < 0.01


def test_check_calibration_is_silent_unless_verbose(capsys):
    surface, _, _, _ = _surface()

    surface.check_calibration(False)
    assert capsys.readouterr().out == ""

    surface.check_calibration(True, tol=0.0)
    out = capsys.readouterr().out
    assert "FWD SWAP RATE" in out
    assert out.count("*") == 15


def test_repr_reports_forward_rates_and_grid():
    surface, _, strike_grid, vol_grid = _surface()

    text = repr(surface)

    assert "SwaptionVolSurface" in text
    assert text.count("FWD SWAP RATE") == 3
    assert text.count("STRIKE, VOL") == strike_grid.size
    assert "STOCK PRICE" not in text
