##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np
import matplotlib.pyplot as plt

import scipy
from scipy.interpolate import splrep

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.global_vars import g_days_in_year
from ...utils.math import scale
from ...utils.helpers import label_to_string

from .curve_fits import CurveFitPolynomial
from .curve_fits import CurveFitNelsonSiegel
from .curve_fits import CurveFitNelsonSiegelSvensson
from .curve_fits import CurveFitBSpline

###############################################################################
# TO DO: CONSTRAIN TAU'S IN NELSON-SIEGEL
###############################################################################


class BondYieldCurve:
    """Class to do fitting of the yield curve and to enable interpolation of
    yields. Because yields assume a flat term structure for each bond, this
    class does not allow discounting to be done and so does not inherit from
    FinDiscountCurve. It should only be used for visualisation and simple
    interpolation but not for full term-structure-consistent pricing."""

    def __init__(
        self, settle_dt: Date, bonds: list, ylds: (np.ndarray, list), curve_fit
    ):
        """Fit the curve to a set of bond yields using the type of curve
        specified. Bounds can be provided if you wish to enforce lower and
        upper limits on the respective model parameters."""

        self.settle_dt = settle_dt
        self.bonds = bonds
        self.ylds = np.array(ylds)
        self.curve_fit = curve_fit

        fit_type = type(self.curve_fit)

        years_to_maturities = []

        for bond in bonds:
            years_to_maturity = (bond.maturity_dt - settle_dt) / g_days_in_year
            years_to_maturities.append(years_to_maturity)

        self.yearsToMaturity = np.array(years_to_maturities)

        if fit_type is CurveFitPolynomial:

            d = curve_fit.power
            coeffs = np.polyfit(self.yearsToMaturity, self.ylds, deg=d)
            curve_fit.coeffs = coeffs

        elif fit_type is CurveFitNelsonSiegel:

            xdata = self.yearsToMaturity
            ydata = self.ylds

            popt, pcov = scipy.optimize.curve_fit(
                curve_fit.interp_yield, xdata, ydata, bounds=curve_fit.bounds
            )

            curve_fit.beta_1 = popt[0]
            curve_fit.beta_2 = popt[1]
            curve_fit.beta_3 = popt[2]
            curve_fit.tau = popt[3]

        elif fit_type is CurveFitNelsonSiegelSvensson:

            xdata = self.yearsToMaturity
            ydata = self.ylds

            popt, pcov = scipy.optimize.curve_fit(
                curve_fit.interp_yield, xdata, ydata, bounds=curve_fit.bounds
            )

            curve_fit.beta_1 = popt[0]
            curve_fit.beta_2 = popt[1]
            curve_fit.beta_3 = popt[2]
            curve_fit.beta_4 = popt[3]
            curve_fit.tau_1 = popt[4]
            curve_fit.tau_2 = popt[5]

        elif fit_type is CurveFitBSpline:

            xdata = self.yearsToMaturity
            ydata = self.ylds

            """ Cubic splines as k=3 """
            spline = splrep(xdata, ydata, k=curve_fit.power, t=curve_fit.knots)
            self.curve_fit.spline = spline

        else:
            raise FinError("Unrecognised curve fit type.")

    ###########################################################################

    def interp_yield(self, maturity_dt: Date):

        if isinstance(maturity_dt, Date):
            t = (maturity_dt - self.settle_dt) / g_days_in_year
        elif isinstance(maturity_dt, list):
            t = maturity_dt
        elif isinstance(maturity_dt, np.ndarray):
            t = maturity_dt
        elif isinstance(maturity_dt, float) or isinstance(
            maturity_dt, np.float64
        ):
            t = maturity_dt
        else:
            raise FinError("Unknown date type.")

        fit = self.curve_fit

        if isinstance(fit, CurveFitPolynomial):
            yld = fit.interp_yield(t)
        elif isinstance(fit, CurveFitNelsonSiegel):
            yld = fit.interp_yield(
                t, fit.beta_1, fit.beta_2, fit.beta_3, fit.tau
            )

        elif isinstance(fit, CurveFitNelsonSiegelSvensson):
            yld = fit.interp_yield(
                t,
                fit.beta_1,
                fit.beta_2,
                fit.beta_3,
                fit.beta_4,
                fit.tau_1,
                fit.tau_2,
            )

        elif isinstance(fit, CurveFitBSpline):
            yld = fit.interp_yield(t)

        return yld

    ###########################################################################

    def plot(self, title, ylabel="Yield To Maturity (%)"):
        """Display yield curve."""

        plt.figure(figsize=(12, 6))
        plt.title(title)
        bond_ylds_scaled = scale(self.ylds, 100.0)
        plt.plot(self.yearsToMaturity, bond_ylds_scaled, "o")
        plt.xlabel("Time to Maturity (years)")
        plt.ylabel(ylabel)

        tmax = np.max(self.yearsToMaturity)
        t = np.linspace(0.0, int(tmax + 0.5), 100)

        yld = self.interp_yield(t)
        yld = scale(yld, 100.0)
        plt.plot(t, yld, label=str(self.curve_fit))
        plt.legend(loc="lower right")
        plt.ylim((min(yld) - 0.3, max(yld) * 1.1))
        plt.grid(True)

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("SETTLEMENT DATE", self.settle_dt)
        s += label_to_string("BOND", self.bonds)
        s += label_to_string("YIELDS", self.ylds)
        s += label_to_string("CURVE FIT", self.curve_fit)
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


###############################################################################
