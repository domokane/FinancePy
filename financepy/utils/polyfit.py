# Fast numba compatible implementation of polyfit
# See https://gist.github.com/kadereub/9eae9cff356bb62cdbd672931e8e5ec4#file-numba_polyfit-py


import numpy as np
import numba

# Goal is to implement a numba compatible polyfit (note does not include
# error handling)

# Define Functions Using Numba
# Idea here is to solve ax = b, using least squares, where a represents
# our coefficients e.g. x**2, x, constants

########################################################################################


@numba.njit("f8[:,:](f8[:], i8)")
def _coeff_mat(x, deg):
    """Construct the coefficient matrix for polynomial fitting."""
    mat = np.zeros(shape=(x.shape[0], deg + 1))
    c = np.ones_like(x)
    mat[:, 0] = c
    mat[:, 1] = x
    if deg > 1:
        for n in range(2, deg + 1):
            mat[:, n] = x**n
    return mat


########################################################################################


@numba.njit("f8[:](f8[:,:], f8[:])")
def _fit_x(a, b):
    """Fit the polynomial coefficients to the data points."""
    # linalg solves ax = b
    det = np.linalg.lstsq(a, b)[0]
    return det


########################################################################################


@numba.njit("f8[:](f8[:], f8[:], i8)")
def fit_poly(x, y, deg):
    """Fit a polynomial of degree `deg` to points (x, y)"""
    a = _coeff_mat(x, deg)
    p = _fit_x(a, y)
    # Reverse order so p[0] is coefficient of highest order
    ret = p[::-1]
    return ret


########################################################################################


@numba.njit(fastmath=True, cache=True)
def eval_polynomial(p, x):
    """
    Compute polynomial p(x) where p is a vector of coefficients, highest
    order coefficient at p[0].  Uses Horner's Method.
    """

    result = np.zeros_like(x)
    for coeff in p:
        result = x * result + coeff
    return result


########################################################################################


# if __name__ == "__main__":

#     # Create Dummy Data and use existing numpy polyfit as test
#     x = np.linspace(0, 2, 20)
#     y = np.cos(x) + 0.3*np.random.rand(20)
#     p = np.poly1d(np.polyfit(x, y, 3))

#     t = np.linspace(0, 2, 200)
#     plt.plot(x, y, 'o', t, p(t), '-')

#     # Now plot using the Numba (amazing) functions
#     p_coeffs = fit_poly(x, y, deg=3)
#     plt.plot(x, y, 'o', t, eval_polynomial(p_coeffs, t), '-')
