###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import copy
import numpy as np
from ...utils.error import FinError
from ...utils.helpers import label_to_string
from ...utils.global_vars import G_SMALL
from scipy.interpolate import BSpline

###############################################################################


class CurveFitMethod:
    """Base class / marker for curve fit methods."""

    def interp_rate(self, t):
        raise NotImplementedError("Subclasses must implement interp_rate().")

    def get_params(self):
        raise NotImplementedError()

    def set_params(self, params):
        raise NotImplementedError()

    def _print(self):
        print(self)

    def clone(self):
        return copy.deepcopy(self)

    @property
    def bounds(self):
        return (-np.inf, np.inf)

###############################################################################
# These are parametric curve fitters that can be used in a number of cases
###############################################################################


class CurveFitPolynomial(CurveFitMethod):
    """Polynomial curve fitting."""

    def __init__(self, power=3):
        self.name = "Polynomial (" + str(power) + ")"
        self.power = power
        n_coeffs = power + 1
        self.coeffs = np.full(n_coeffs, 0.03)
        self._bounds = (-np.inf, np.inf)

    def interp_rate(self, t):
        # I store coefficients with lowest power first but numpy wants
        # the highest power first so I need to reverse the order
        coeffs = self.coeffs[::-1]
        t = np.asarray(t, dtype=float)
        yld = np.polyval(coeffs, t)
        return yld

    def get_params(self):
        return np.array(self.coeffs, dtype=float)

    def set_params(self, params):
        self.coeffs = np.array(params, dtype=float)

    @property
    def bounds(self):
        return self._bounds

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("Power", self.power)

        for c in self.coeffs:
            s += label_to_string("Coefficient", c)

        return s

###############################################################################

class CurveFitNelsonSiegel(CurveFitMethod):
    """Nelson-Siegel parametric fit."""

    def __init__(self, beta1=0.03, beta2=-0.02, beta3=0.02, tau=2.0, bounds=None):
        self.name = "Nelson-Siegel"
        self.beta_1 = beta1
        self.beta_2 = beta2
        self.beta_3 = beta3
        self.tau = tau

        # Fairly permissive bounds. Only tau is restricted to 0.5-100.
        if bounds is None:
            bounds = [(-5.0, -5.0, -5.0, 0.01),
                      (10.0, 10.0, 10.0, 100.0)]

        self._bounds = bounds

    def get_params(self):
        return np.array([self.beta_1, self.beta_2, self.beta_3, self.tau], dtype=float)

    def set_params(self, params):
        self.beta_1, self.beta_2, self.beta_3, self.tau = params

    @property
    def bounds(self):
        return self._bounds

    def interp_rate(self, t, beta_1=None, beta_2=None, beta_3=None, tau=None):
        # This rate is not specific.

        t = np.asarray(t, dtype=float)
        t = np.maximum(t, 1e-10)

        if beta_1 is None: beta_1 = self.beta_1
        if beta_2 is None: beta_2 = self.beta_2
        if beta_3 is None: beta_3 = self.beta_3
        if tau is None:    tau = self.tau

        if tau <= G_SMALL:
            raise FinError("tau must be positive")

        theta = t / tau
        exp_term = np.exp(-theta)
        yld = beta_1
        yld += beta_2 * (1.0 - exp_term) / theta
        yld += beta_3 * ((1.0 - exp_term) / theta - exp_term)
        return yld

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("beta_1", self.beta_1)
        s += label_to_string("beta_2", self.beta_2)
        s += label_to_string("beta_3", self.beta_3)
        s += label_to_string("tau", self.tau)
        return s

###############################################################################

class CurveFitSvensson(CurveFitMethod):
    """Svensson (extended Nelson-Siegel) parametric fit."""

    def __init__(
        self,
        beta1=0.03,
        beta2=-0.02,
        beta3=0.02,
        beta4=0.01,
        tau1=2.0,
        tau2=5.0,
        bounds=None,
        ):

        self.name = "Svensson"
        self.beta_1 = beta1
        self.beta_2 = beta2
        self.beta_3 = beta3
        self.beta_4 = beta4
        self.tau_1 = tau1
        self.tau_2 = tau2

        if bounds is None:
            bounds = [(-5.0, -5.0, -5.0, -5.0, 0.01, 0.1),
                      (10.0, 10.0, 10.0, 10.0, 10.0, 100.0)]
        self._bounds = bounds

    def get_params(self):
        return np.array(
            [self.beta_1,
             self.beta_2,
             self.beta_3,
             self.beta_4,
             self.tau_1,
             self.tau_2],
            dtype=float,
        )

    def set_params(self, params):
        (
            self.beta_1,
            self.beta_2,
            self.beta_3,
            self.beta_4,
            self.tau_1,
            self.tau_2,
        ) = params

    @property
    def bounds(self):
        return self._bounds

    def interp_rate(
        self,
        t,
        beta_1=None,
        beta_2=None,
        beta_3=None,
        beta_4=None,
        tau_1=None,
        tau_2=None,
    ):
        # This rate is not specific. In some cases will be the bond
        # yield to maturity. In other cases it will be the zero rate
        # or forward rate

        t = np.asarray(t, dtype=float)
        t = np.maximum(t, 1e-10)

        if beta_1 is None: beta_1 = self.beta_1
        if beta_2 is None: beta_2 = self.beta_2
        if beta_3 is None: beta_3 = self.beta_3
        if beta_4 is None: beta_4 = self.beta_4
        if tau_1 is None: tau_1 = self.tau_1
        if tau_2 is None: tau_2 = self.tau_2

        if tau_1 <= G_SMALL or tau_2 <= G_SMALL:
            raise FinError("tau1 and tau2 must be positive")

        theta1 = t / tau_1
        theta2 = t / tau_2
        exp_term1 = np.exp(-theta1)
        exp_term2 = np.exp(-theta2)
        yld = beta_1
        yld += beta_2 * (1.0 - exp_term1) / theta1
        yld += beta_3 * ((1.0 - exp_term1) / theta1 - exp_term1)
        yld += beta_4 * ((1.0 - exp_term2) / theta2 - exp_term2)
        return yld

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("beta_1", self.beta_1)
        s += label_to_string("beta_2", self.beta_2)
        s += label_to_string("beta_3", self.beta_3)
        s += label_to_string("beta_4", self.beta_4)
        s += label_to_string("tau_1", self.tau_1)
        s += label_to_string("tau_2", self.tau_2)
        return s


###############################################################################

class CurveFitBSpline(CurveFitMethod):
    """B-Spline curve fitting."""

    def __init__(self, power=3, knot_years=None, t_max=50.0):
        self.name = "B-Spline"
        self.power = power
        self.t_min = 0.0
        self.t_max = t_max

        if knot_years is None:
            knot_years = [1.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0]

        self.knot_years = np.asarray(knot_years, dtype=float)

        if np.any(self.knot_years <= 0.0):
            raise FinError("Knot years must be positive.")

        if np.any(self.knot_years >= t_max):
            raise FinError(f"Knot years must be less than t_max={t_max}.")

        # Normalise knots to [0, 1]
        self.knots = self.knot_years / t_max

        k = power

        self.t = np.concatenate((
            np.full(k + 1, 0.0),
            self.knots,
            np.full(k + 1, 1.0),
        ))

        n_coeffs = len(self.t) - k - 1
        self.coeffs = np.full(n_coeffs, 0.03)
        self._spline = BSpline(self.t, self.coeffs, self.power, extrapolate=True)
        self._bounds = (-np.inf, np.inf)

    def get_params(self):
        return np.array(self.coeffs, dtype=float)

    def set_params(self, params):
        self.coeffs = np.array(params, dtype=float)
        self._spline = BSpline(self.t, self.coeffs, self.power, extrapolate=True)

    @property
    def bounds(self):
        return self._bounds

    def interp_rate(self, t):
        t = np.clip(np.asarray(t, dtype=float), self.t_min, self.t_max)
        return self._spline(t)

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("Power", self.power)
        s += label_to_string("t_min", self.t_min)
        s += label_to_string("t_max", self.t_max)
        s += label_to_string("Knots", self.knots)
        for i, c in enumerate(self.coeffs):
            s += label_to_string(f"Coefficient {i}", c)
        return s

###############################################################################
