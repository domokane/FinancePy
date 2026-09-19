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

from financepy.models.gbm_process_simulator import get_assets_paths_times
from financepy.utils.math import corr_matrix_generator


########################################################################################


def test_fin_gbm_process():

    num_assets = 3
    num_paths = 6
    num_time_steps = 1
    t = 1.0
    mus = 0.03 * np.ones(num_assets)
    stock_prices = 100.0 * np.ones(num_assets)
    volatilities = 0.2 * np.ones(num_assets)
    rho = 0.8
    corr_matrix = corr_matrix_generator(rho, num_assets)
    seed = 1912

    times, paths = get_assets_paths_times(
        num_assets,
        num_paths,
        num_time_steps,
        t,
        mus,
        stock_prices,
        volatilities,
        corr_matrix,
        seed,
    )


########################################################################################

test_fin_gbm_process()
