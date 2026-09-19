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

from financepy.models.vasicek_mc import zero_price, zero_price_mc



########################################################################################


def test_fin_model_rates_vasicek():

    r0 = 0.05
    a = 0.10
    b = 0.05
    sigma = 0.05
    t = 5.0

    p = zero_price(r0, a, b, sigma, t)

    num_paths = 1000
    dt = 0.02
    seed = 1968

    print("TIME", "T", "P", "P_MC", "P_MC2")

    for t in np.linspace(0, 10, 21):
        start = time.time()
        p_mc = zero_price_mc(r0, a, b, sigma, t, dt, num_paths, seed)
        p_mc2 = zero_price_mc(r0, a, b, sigma, t, dt, 10 * num_paths, seed)
        p = zero_price(r0, a, b, sigma, t)
        end = time.time()
        elapsed = end - start
        print(elapsed, t, p, p_mc, p_mc2)


########################################################################################

test_fin_model_rates_vasicek()
