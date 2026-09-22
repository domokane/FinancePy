##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import matplotlib.pyplot as plt

from typing import Union

import numpy as np

import scipy
from scipy.optimize import least_squares

from ...utils.format_graphs import *
from ...utils.error import FinError
from ...utils.date import Date
from ...utils.day_count import DayCountTypes
from ...utils.math import scale
from ...utils.helpers import label_to_string
from ...utils.helpers import times_from_dates

from .curve_fits import CurveFitTypes
from .curve_fits import CurveFitCubicPolynomial
from .curve_fits import CurveFitQuarticPolynomial
from .curve_fits import CurveFitQuinticPolynomial
from .curve_fits import CurveFitNelsonSiegel
from .curve_fits import CurveFitNelsonSiegelSvensson
from .curve_fits import CurveFitBSpline


class BondParametricYieldCurve:
    """Fit and interpolate a bond yield curve."""

    def __init__(
        self,
        settle_dt: Date,
        bonds: list,
        ylds: Union[np.ndarray, list],
        curve_fit_type: CurveFitTypes,
        curve_dc_type: DayCountTypes = DayCountTypes.ACT_365F,
    ):

        self.settle_dt = settle_dt
        self.bonds = bonds
        self.ylds = np.asarray(
            ylds,
            dtype=float,
        )

        if not isinstance(curve_fit_type, CurveFitTypes):
            raise FinError("Invalid curve fit type.")

        self.curve_fit_type = curve_fit_type

        ######################################################################
        # Create curve fitter
        ######################################################################

        if curve_fit_type == CurveFitTypes.CUBIC_POLYNOMIAL:
            self.curve_fit = CurveFitCubicPolynomial()

        elif curve_fit_type == CurveFitTypes.QUARTIC_POLYNOMIAL:
            self.curve_fit = CurveFitQuarticPolynomial()

        elif curve_fit_type == CurveFitTypes.QUINTIC_POLYNOMIAL:
            self.curve_fit = CurveFitQuinticPolynomial()

        elif curve_fit_type == CurveFitTypes.NELSON_SIEGEL:
            self.curve_fit = CurveFitNelsonSiegel()

        elif curve_fit_type == CurveFitTypes.NELSON_SIEGEL_SVENSSON:
            self.curve_fit = CurveFitNelsonSiegelSvensson()

        elif curve_fit_type == CurveFitTypes.BSPLINE:
            self.curve_fit = CurveFitBSpline()

        else:
            raise FinError("Unrecognised curve fit type.")

        ######################################################################

        if not isinstance(curve_dc_type, DayCountTypes):
            raise FinError("Invalid time day count type.")

        self.curve_dc_type = curve_dc_type

        ######################################################################
        # Calculate times to maturity
        ######################################################################

        years_to_maturities = []

        for bond in bonds:

            t = times_from_dates(
                settle_dt,
                bond.maturity_dt,
                self.curve_dc_type,
            )

            years_to_maturities.append(t)

        self.years_to_maturity = np.asarray(
            years_to_maturities,
            dtype=float,
        )

        self.t_max = max(
            np.max(self.years_to_maturity),
            1.0e-8,
        )

        tdata = self.years_to_maturity
        ylds = self.ylds

        ######################################################################
        # Polynomial fit
        ######################################################################

        if curve_fit_type in (
            CurveFitTypes.CUBIC_POLYNOMIAL,
            CurveFitTypes.QUARTIC_POLYNOMIAL,
            CurveFitTypes.QUINTIC_POLYNOMIAL,
        ):

            self.curve_fit.t_scale = self.t_max

            xdata = tdata / self.curve_fit.t_scale

            degree = self.curve_fit.power

            coeffs_high_first = np.polyfit(
                xdata,
                ylds,
                deg=degree,
            )

            self.curve_fit.coeffs = coeffs_high_first[::-1]

        ######################################################################
        # Nelson-Siegel fit
        ######################################################################

        elif curve_fit_type == CurveFitTypes.NELSON_SIEGEL:

            popt, _ = scipy.optimize.curve_fit(
                self.curve_fit.interp_rate,
                tdata,
                ylds,
                p0=self.curve_fit.get_params(),
                bounds=self.curve_fit.bounds,
                maxfev=10000,
            )

            self.curve_fit.set_params(popt)

        ######################################################################
        # Svensson fit
        ######################################################################

        elif curve_fit_type == CurveFitTypes.NELSON_SIEGEL_SVENSSON:

            popt, _ = scipy.optimize.curve_fit(
                self.curve_fit.interp_rate,
                tdata,
                ylds,
                p0=self.curve_fit.get_params(),
                bounds=self.curve_fit.bounds,
                maxfev=10000,
            )

            self.curve_fit.set_params(popt)

        ######################################################################
        # B-Spline fit
        ######################################################################

        elif curve_fit_type == CurveFitTypes.BSPLINE:

            def residuals(params):

                self.curve_fit.set_params(params)

                return self.curve_fit.interp_rate(tdata) - ylds

            result = least_squares(
                residuals,
                self.curve_fit.get_params(),
                bounds=self.curve_fit.bounds,
                xtol=1.0e-10,
                ftol=1.0e-10,
                gtol=1.0e-10,
                max_nfev=1000,
            )

            if not result.success:
                raise FinError(result.message)

            self.curve_fit.set_params(result.x)

        else:

            raise FinError("Unrecognised curve fit type.")

    ##########################################################################

    def interp_yield(self, maturity_dt):
        """Interpolate yield."""

        if isinstance(maturity_dt, Date):

            t = times_from_dates(
                self.settle_dt,
                maturity_dt,
                self.curve_dc_type,
            )

        elif isinstance(
            maturity_dt,
            (
                list,
                np.ndarray,
                float,
                int,
                np.floating,
                np.integer,
            ),
        ):

            t = maturity_dt

        else:

            raise FinError("Unknown date type.")

        return self.curve_fit.interp_rate(t)

    ##########################################################################

    def errors(self):
        """Return RMS and maximum fit errors in basis points."""

        ylds = self.ylds
        times = self.years_to_maturity

        y_fit = self.curve_fit.interp_rate(times)

        res = ylds - y_fit

        mean_err = np.sqrt(np.mean(res * res))

        max_err = np.max(np.abs(res))

        bp = 10000.0

        return (
            mean_err * bp,
            max_err * bp,
        )

    ##########################################################################

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

        plt.figure()

        plot_title = title + " - " + self.curve_fit.name

        plt.title(plot_title)

        if times is None:

            tmax = np.max(self.years_to_maturity)

            n_points = max(
                int(12 * tmax),
                2,
            )

            times = np.linspace(
                0.0,
                tmax,
                n_points,
            )

        else:

            times = np.asarray(
                times,
                dtype=float,
            )

        if np.any(times < 0.0):

            raise FinError("Plot times must be non-negative.")

        times = np.maximum(
            times,
            1.0e-8,
        )

        bond_ylds_scaled = scale(
            self.ylds,
            100.0,
        )

        plt.plot(
            self.years_to_maturity,
            bond_ylds_scaled,
            "o",
            label="Bond Yields",
        )

        ytm = self.interp_yield(times)

        plt.xlabel("Time to Maturity (years)")

        plt.ylabel("Yield (%)")

        plt.plot(
            times,
            ytm * 100.0,
            label=self.curve_fit.name,
        )

        plt.legend(loc="lower right")

        if ymin is not None and ymax is not None:

            plt.ylim(
                ymin,
                ymax,
            )

        plt.xlim(
            np.min(times),
            np.max(times),
        )

        plt.grid(
            True,
            alpha=0.3,
        )

        if filename is not None:

            plt.savefig(
                filename,
                bbox_inches="tight",
                pad_inches=0.02,
            )

        plt.show()
        plt.close()

    ##########################################################################

    def __repr__(self):

        s = label_to_string(
            "OBJECT TYPE",
            type(self).__name__,
        )

        s += label_to_string(
            "SETTLEMENT DATE",
            self.settle_dt,
        )

        s += label_to_string(
            "YIELDS",
            self.ylds,
        )

        s += label_to_string(
            "CURVE FIT TYPE",
            self.curve_fit_type,
        )

        s += label_to_string(
            "CURVE FIT",
            self.curve_fit,
        )

        return s

    ##########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""

        print(self)


###############################################################################
