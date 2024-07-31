import pytest
import numpy as np

from financepy.utils.global_vars import *
from financepy.utils.tension_spline import TensionSpline
from financepy.utils.global_vars import g_small



@pytest.mark.parametrize('n_knots', [1, 2, 3, 4, 5])
def test_TensorSpline_Recovers_Inputs(n_knots):
    x = np.arange(n_knots)
    y = x*0.1 + np.arange(n_knots) % 2
    sigma = 1.0
    ts = TensionSpline(x, y, sigma)
    y_out = ts(x)
    relerr = np.linalg.norm(y-y_out)/(np.linalg.norm(y) + g_small)
    assert relerr <= g_small


def test_TensorSpline_Values():
    n_knots = 4
    x = np.arange(n_knots)
    y = x*0.1 + np.arange(n_knots) % 2
    sigma = 1.0
    ts = TensionSpline(x, y, sigma)

    xs = np.linspace(x[0]-1, x[-1]+1, 2*n_knots+1, endpoint=True)
    y_out = ts(xs)

    y_expected = [0.0, 0.0, 0.42195359, 1.0955618, 0.65, 0.2044382, 0.87804641, 1.3, 1.3, ]
    relerr = np.linalg.norm(y_expected-y_out)/np.linalg.norm(y_expected)
    assert relerr <= g_small*1e4  # some rounding errors expected


# if __name__ == '__main__':
    # test_TensorSpline_Recovers_Inputs(2)
    # test_TensorSpline_Values()
