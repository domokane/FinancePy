import numpy as np
from scipy.linalg import solve_banded


class TensionSpline(object):
    """
    Implement exponential tension spline following [AP] Sec 6.A.3

    [AP] Andersen, Piterbarg. Interest Rate Modeling, 2010
    """

    def __init__(self, x, y, sigma):
        self._x = np.atleast_1d(x).astype(float)
        self._y = np.atleast_1d(y).astype(float)
        self._sigma = max(sigma, 1e-2)
        self.calculate_coefs()

    def validate_inputs(self):
        assert len(self._x) == len(
            self._y
        ), "x and y should be the same length"
        assert len(self._x) > 0, "Empty arrays are not allowed"
        assert len(self._x) == 1 or np.min(
            np.diff(self._x) > 1e-8
        ), "x should be strictly increasing"

    def calculate_coefs(self):

        self.validate_inputs()

        N = len(self._x)
        sig = self._sigma

        # prelim calcs
        # TODO: cache 1/h, 1/sh etc not h, sh for efficiency
        dy = np.diff(self._y)
        self._h = np.diff(self._x)
        self._sh = np.sinh(sig * self._h)
        self._ch = np.cosh(sig * self._h)

        # useful views
        hl = self._h[:-1]
        hr = self._h[1:]
        shl = self._sh[:-1]
        shr = self._sh[1:]
        chl = self._ch[:-1]
        chr = self._ch[1:]

        # set up a tri-diagona matrix in the form needed for solve_baded
        ab = np.zeros((3, N))

        # main diagonal
        ab[1, 0] = 1
        ab[1, -1] = 1
        ab[1, 1:-1] = sig * (chl / shl + chr / shr) - 1 / hl - 1 / hr

        # lower diagonal
        ab[2, :-2] = 1 / hl - sig / shl

        # upper diagonal
        ab[0, 2:] = 1 / hr - sig / shr

        # rhs
        b = np.zeros(N)
        b[1:-1] = (dy[1:] / hr - dy[:-1] / hl) * sig * sig

        self._ypp = solve_banded((1, 1), ab, b)

    def __call__(self, xs):
        xs = np.atleast_1d(xs).astype(float)

        # rename for brevity
        h = self._h
        sh = self._sh
        sig = self._sigma
        sig2 = sig * sig

        ids = np.searchsorted(self._x, xs, side="left")
        out = np.zeros_like(xs)
        for i, x in enumerate(xs):
            idx = ids[i]

            # handle extrapolation first
            if idx == 0:
                out[i] = self._y[0]
                continue
            if idx == len(self._x):
                out[i] = self._y[-1]
                continue

            # main calc. self._x[idx-1] < x <= self._x[idx]
            xl = self._x[idx - 1]
            xr = self._x[idx]
            v1 = (
                (np.sinh(sig * (xr - x)) / sh[idx - 1] - (xr - x) / h[idx - 1])
                * self._ypp[idx - 1]
                / sig2
            )
            v2 = (
                (np.sinh(sig * (x - xl)) / sh[idx - 1] - (x - xl) / h[idx - 1])
                * self._ypp[idx]
                / sig2
            )
            v3 = self._y[idx - 1] * (xr - x) / h[idx - 1]
            v4 = self._y[idx] * (x - xl) / h[idx - 1]
            out[i] = v1 + v2 + v3 + v4

        return out
