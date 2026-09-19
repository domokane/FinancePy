# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path
_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause
_install_double_click_pause()

# import numpy as np



# def test_FinPieceCurve():
#    times = np.linspace(0.0, 1.0, 5)
#    values = np.ones(5) * 0.05
#    flat_curve = FinPieceCurve(times,values)
#    dfs = flat_curve.df(times, 0)
#    print(dfs)
#    dfs = flat_curve.df(times, 1)
#    print(dfs)
#    dfs = flat_curve.df(times, 2)
#    print(dfs)
#    dfs = flat_curve.df(times, -1)
#    print(dfs)


# test_FinPieceCurve()
