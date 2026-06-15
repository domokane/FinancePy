##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import matplotlib.pyplot as plt

from typing import Union

import numpy as np

import scipy
from scipy.interpolate import splrep

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.day_count import DayCountTypes
from ...utils.math import scale
from ...utils.helpers import label_to_string
from ...utils.helpers import times_from_dates

from .curve_fits import CurveFitPolynomial
from .curve_fits import CurveFitNelsonSiegel
from .curve_fits import CurveFitNelsonSiegelSvensson
from .curve_fits import CurveFitBSpline


class BondFittedYieldCurve:
    """Class to do fitting of the yield curve and to enable interpolation of
    yields. Because yields assume a flat term structure for each bond, this
    class does not allow discounting to be done and so does not inherit from
    DiscountCurve. It should only be used for visualisation and simple
    interpolation but not for full term-structure-consistent pricing."""

    def __init__(
        self,
        settle_dt: Date,
        bonds: list,
        ylds: Union[np.ndarray, list],
        curve_fit,
        time_dc_type: DayCountTypes = DayCountTypes.ACT_365F,
    ):
        """Fit the curve to a set of bond yields using the type of curve
        specified. Bounds can be provided if you wish to enforce lower and
        upper limits on the respective model parameters."""

        print("Class BondFittedYieldCurve is Deprecated.")
        print("Use BondParametricYieldCurve in market/curves")

        self.settle_dt = settle_dt
        self.bonds = bonds
        self.ylds = np.array(ylds)
        self.curve_fit = curve_fit

        fit_type = type(self.curve_fit)

        if not isinstance(time_dc_type, DayCountTypes):
            raise FinError("Invalid time day count type.")

        self.time_dc_type = time_dc_type

        years_to_maturities = []

        for bond in bonds:

            years_to_maturity = times_from_dates(
                settle_dt,
                bond.maturity_dt,
                self.time_dc_type,
            )
            years_to_maturities.append(years_to_maturity)

        self.years_to_maturity = np.array(years_to_maturities)

        if fit_type is CurveFitPolynomial:

            d = curve_fit.power
            # Highest powers are first
            coeffs1 = np.polyfit(self.years_to_maturity, self.ylds, deg=d)
            # Lowest powers are first
            coeffs2 = coeffs1[::-1]
            curve_fit.coeffs = coeffs2

        elif fit_type is CurveFitNelsonSiegel:

            xdata = self.years_to_maturity
            ydata = self.ylds

            popt, _ = scipy.optimize.curve_fit(
                curve_fit.interp_rate, xdata, ydata, bounds=curve_fit.bounds
            )

            curve_fit.beta_1 = popt[0]
            curve_fit.beta_2 = popt[1]
            curve_fit.beta_3 = popt[2]
            curve_fit.tau = popt[3]

        elif fit_type is CurveFitNelsonSiegelSvensson:

            xdata = self.years_to_maturity
            ydata = self.ylds

            popt, _ = scipy.optimize.curve_fit(
                curve_fit.interp_rate, xdata, ydata, bounds=curve_fit.bounds
            )

            curve_fit.beta_1 = popt[0]
            curve_fit.beta_2 = popt[1]
            curve_fit.beta_3 = popt[2]
            curve_fit.beta_4 = popt[3]
            curve_fit.tau_1 = popt[4]
            curve_fit.tau_2 = popt[5]

        elif fit_type is CurveFitBSpline:

            xdata = self.years_to_maturity
            ydata = self.ylds

            """ Cubic splines as k=3 """
            spline = splrep(xdata, ydata, k=curve_fit.power, t=curve_fit.knots)
            self.curve_fit.spline = spline

        else:
            raise FinError("Unrecognised curve fit type.")

    ###########################################################################

    def interp_yield(self, maturity_dt: Date):
        """Interpolate yield"""
        if isinstance(maturity_dt, Date):
            t = times_from_dates(self.settle_dt, maturity_dt, self.time_dc_type)
        elif isinstance(maturity_dt, list):
            t = maturity_dt
        elif isinstance(maturity_dt, np.ndarray):
            t = maturity_dt
        elif isinstance(maturity_dt, float) or isinstance(maturity_dt, np.float64):
            t = maturity_dt
        else:
            raise FinError("Unknown date type.")

        fit = self.curve_fit
        yld = None

        if isinstance(fit, CurveFitPolynomial):
            yld = fit.interp_rate(t)
        elif isinstance(fit, CurveFitNelsonSiegel):
            yld = fit.interp_rate(t, fit.beta_1, fit.beta_2, fit.beta_3, fit.tau)

        elif isinstance(fit, CurveFitNelsonSiegelSvensson):
            yld = fit.interp_rate(
                t,
                fit.beta_1,
                fit.beta_2,
                fit.beta_3,
                fit.beta_4,
                fit.tau_1,
                fit.tau_2,
            )

        elif isinstance(fit, CurveFitBSpline):
            yld = fit.interp_rate(t)

        return yld

    ###########################################################################

    def plot(
        self,
        title,
        times: np.ndarray = None,
        ymin: float = None,
        ymax: float = None,
        filename: str = None,
    ):
        """Display yield curve."""

        plt.rcParams.update(
            {
                "lines.linewidth": 3,
                "font.size": 14,
                "axes.labelsize": 14,
                "axes.titlesize": 16,
                "legend.fontsize": 14,
            }
        )

        plt.figure(figsize=(12, 6))

        title = title + " - " + self.curve_fit.name
        plt.title(title)

        if times is None:
            tmax = np.max(self.years_to_maturity)
            times = np.linspace(0.0, tmax, int(12 * tmax))
        else:
            times = np.asarray(times, dtype=float)

        if np.any(times < 0.0):
            raise FinError("Plot times must be strictly positive.")

        times = np.maximum(times, 1e-8)

        bond_ylds_scaled = scale(self.ylds, 100.0)

        # Plot actual bond yields
        plt.plot(self.years_to_maturity, bond_ylds_scaled, "o", label="Bond Yields")

        ytm = self.interp_yield(times)

        plt.xlabel("Time to Maturity (years)")
        plt.ylabel("Yield (%)")

        plt.plot(times, ytm * 100, label=str(self.curve_fit.name))
        plt.legend(loc="lower right")

        if ymin is not None and ymax is not None:
            plt.ylim(ymin, ymax)

        plt.xlim(np.min(times), np.max(times))
        plt.grid(True, alpha=0.3)
        plt.tight_layout()

        if filename is not None:
            plt.savefig(filename, bbox_inches="tight", pad_inches=0.02)

        plt.show()

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
