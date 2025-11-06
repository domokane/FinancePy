##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from typing import Union

import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import splrep

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.global_vars import G_DAYS_IN_YEAR
from ...utils.math import scale
from ...utils.helpers import label_to_string

from .curve_fits import CurveFitPolynomial
from .curve_fits import CurveFitNelsonSiegel
from .curve_fits import CurveFitNelsonSiegelSvensson
from .curve_fits import CurveFitBSpline

########################################################################################
# TO DO: CONSTRAIN TAU'S IN NELSON-SIEGEL
########################################################################################


class BondYieldCurve:
    """Class to do fitting of the yield curve and to enable interpolation of
    yields. Because yields assume a flat term structure for each bond, this
    class does not allow discounting to be done and so does not inherit from
    DiscountCurve. It should only be used for visualisation and simple
    interpolation but not for full term-structure-consistent pricing."""

    def __init__(
        self,
        settlement_date: Date,
        bonds: list,
        ylds: Union[np.ndarray, list],
        curve_fit,
    ):
        """Fit the curve to a set of bond yields using the type of curve
        specified. Bounds can be provided if you wish to enforce lower and
        upper limits on the respective model parameters."""

        self._settlement_date = settlement_date
        self._bonds = bonds
        self._ylds = np.array(ylds)
        self._curve_fit = curve_fit

        fit_type = type(self._curve_fit)
        fit = self._curve_fit

        years_to_maturities = []
        for bond in bonds:
            years_to_maturity = (bond._maturity_dt - settlement_date) / G_DAYS_IN_YEAR
            years_to_maturities.append(years_to_maturity)
        self._years_to_maturity = np.array(years_to_maturities)

        if fit_type is CurveFitPolynomial:

            d = fit._power
            coeffs = np.polyfit(self._years_to_maturity, self._ylds, deg=d)
            fit._coeffs = coeffs

        elif fit_type is CurveFitNelsonSiegel:

            xdata = self._years_to_maturity
            ydata = self._ylds

            popt, pcov = curve_fit(
                self._curve_fit._interpolated_yield,
                xdata,
                ydata,
                bounds=fit._bounds,
            )

            fit._beta1 = popt[0]
            fit._beta2 = popt[1]
            fit._beta3 = popt[2]
            fit._tau = popt[3]

        elif fit_type is CurveFitNelsonSiegelSvensson:

            xdata = self._years_to_maturity
            ydata = self._ylds

            popt, pcov = curve_fit(
                self._curve_fit._interpolated_yield,
                xdata,
                ydata,
                bounds=fit._bounds,
            )

            fit._beta1 = popt[0]
            fit._beta2 = popt[1]
            fit._beta3 = popt[2]
            fit._beta4 = popt[3]
            fit._tau1 = popt[4]
            fit._tau2 = popt[5]

        elif fit_type is CurveFitBSpline:

            xdata = self._years_to_maturity
            ydata = self._ylds

            """ Cubic splines as k=3 """
            spline = splrep(xdata, ydata, k=fit._power, t=fit._knots)
            fit._spline = spline

        else:
            raise FinError("Unrecognised curve fit type.")

    ####################################################################################

    @property
    def curve_fit(self):
        """accessor function for curve_fit"""
        return self._curve_fit

    ####################################################################################

    def interpolated_yield(self, maturity_dt: Date):
        """Interpolates the yield for a given maturity date."""
        if isinstance(maturity_dt, Date):
            t = (maturity_dt - self._settlement_date) / G_DAYS_IN_YEAR
        elif isinstance(maturity_dt, list):
            t = maturity_dt
        elif isinstance(maturity_dt, np.ndarray):
            t = maturity_dt
        elif isinstance(maturity_dt, (float, np.float64)):
            t = maturity_dt
        else:
            raise FinError("Unknown date type.")

        fit = self._curve_fit
        yld = None

        if isinstance(fit, CurveFitPolynomial):
            yld = fit._interpolated_yield(t)
        elif isinstance(fit, CurveFitNelsonSiegel):
            yld = fit._interpolated_yield(
                t, fit._beta1, fit._beta2, fit._beta3, fit._tau
            )

        elif isinstance(fit, CurveFitNelsonSiegelSvensson):
            yld = fit._interpolated_yield(
                t,
                fit._beta1,
                fit._beta2,
                fit._beta3,
                fit._beta4,
                fit._tau1,
                fit._tau2,
            )

        elif isinstance(fit, CurveFitBSpline):
            yld = fit._interpolated_yield(t)

        return yld

    ####################################################################################

    def plot(self, title, ylabel="Yield To Maturity (%)"):
        """Display yield curve."""

        plt.figure(figsize=(12, 6))
        plt.title(title)
        bond_ylds_scaled = scale(self._ylds, 100.0)
        plt.plot(self._years_to_maturity, bond_ylds_scaled, "o")
        plt.xlabel("Time to Maturity (years)")
        plt.ylabel(ylabel)

        tmax = np.max(self._years_to_maturity)
        t = np.linspace(0.0, int(tmax + 0.5), 100)

        yld = self.interpolated_yield(t)
        yld = scale(yld, 100.0)
        plt.plot(t, yld, label=str(self._curve_fit))
        plt.legend(loc="lower right")
        plt.ylim((min(yld) - 0.3, max(yld) * 1.1))
        plt.grid(True)

    ####################################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("SETTLEMENT DATE", self._settlement_date)
        s += label_to_string("BOND", self._bonds)
        s += label_to_string("YIELDS", self._ylds)
        s += label_to_string("CURVE FIT", self._curve_fit)
        return s

    ####################################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


##############################################################################
