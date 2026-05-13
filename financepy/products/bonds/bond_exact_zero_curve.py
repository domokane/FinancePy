##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import matplotlib.pyplot as plt

import numpy as np
from scipy import optimize

from ...utils.date import Date
from ...utils.math import scale, test_monotonicity
from ...utils.global_vars import G_DAYS_IN_YEAR
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.helpers import input_time
from ...utils.helpers import table_to_string
from ...market.curves.interpolator import InterpTypes, interpolate
from ...utils.error import FinError
from ...utils.frequency import annual_frequency, FrequencyTypes
from ...market.curves.discount_curve import DiscountCurve
from ...utils.helpers import label_to_string

########################################################################################


def _f(df, *args):
    curve = args[0]
    value_dt = args[1]
    bond = args[2]
    mkt_clean_price = args[3]

    curve.set_last_df(df)
    bond_discount_price = bond.clean_price_from_discount_curve(value_dt, curve)
    obj_fn = bond_discount_price - mkt_clean_price

    return obj_fn


########################################################################################


class BondExactZeroCurve(DiscountCurve):
    """Class to do bootstrap exact fitting of the bond zero rate curve."""

    def __init__(
        self,
        value_dt: Date,
        bonds: list,
        clean_prices: list,
        interp_type: InterpTypes = InterpTypes.FLAT_FWD_RATES,
        time_dc_type: DayCountTypes = DayCountTypes.ACT_365F,
    ):
        """Fit a discount curve to a set of bond prices using the type of
        curve specified."""

        if len(bonds) != len(clean_prices):
            raise FinError("Num bonds does not equal number of prices.")

        self.time_dc_type = time_dc_type
        self.settle_dt = value_dt
        self.value_dt = value_dt
        self.bonds = bonds
        self.clean_prices = np.array(clean_prices, dtype=float)
        self._interp_type = interp_type

        times = []
        for bond in self.bonds:
            t_mat = (bond.maturity_dt - self.settle_dt) / G_DAYS_IN_YEAR
            times.append(t_mat)

        times = np.array(times, dtype=float)
        if test_monotonicity(times) is False:
            raise FinError("Times are not sorted in increasing order")

        self.years_to_maturity = np.array(times)
        self._bootstrap_zero_rates()

    ###########################################################################

    def _bootstrap_zero_rates(self):

        self._times = np.array([0.0], dtype=float)
        self._dfs = np.array([1.0], dtype=float)

        for i in range(0, len(self.bonds)):

            bond = self.bonds[i]
            maturity_dt = bond.maturity_dt
            clean_price = self.clean_prices[i]
            t_mat = (maturity_dt - self.settle_dt) / G_DAYS_IN_YEAR

            # Let's give it a good starting guess
            df = np.exp(-t_mat * bond.cpn)

            argtuple = (self, self.settle_dt, bond, clean_price)
            self._times = np.append(self._times, t_mat)
            self._dfs = np.append(self._dfs, df)

            df_star = optimize.newton(
                _f,
                x0=df,
                fprime=None,
                args=argtuple,
                tol=1e-8,
                maxiter=100,
                fprime2=None,
            )

            if i > 0 and df_star > self._dfs[i]:
                raise FinError(
                    "Bootstrapped discount factors must be non-increasing"
                )

            if df_star <= 0.0:
                raise FinError("Bootstrapped discount factor must be positive")

            self.set_last_df(df_star)

    ###########################################################################

    def df_t(self, t):
        """Discount factor from time or vector of times."""
        return interpolate(t, self._times, self._dfs, self._interp_type.value)

    # ###########################################################################

    # def zero_rate(
    #     self,
    #     dt: Date,
    #     freq_type: FrequencyTypes = FrequencyTypes.CONTINUOUS,
    # ):
    #     """Calculate the zero rate to maturity date."""
    #     t = input_time(dt, self)
    #     f = annual_frequency(freq_type)
    #     df = self.df_t(t)

    #     if f == 0:  # Simple interest
    #         zero_rate = (1.0 / df - 1.0) / t
    #     elif f == -1:  # Continuous
    #         zero_rate = -np.log(df) / t
    #     else:
    #         zero_rate = (df ** (-1.0 / (f * t)) - 1) * f
    #     return zero_rate

    # ###########################################################################

    # def fwd(self, dt):
    #     """Calculate the continuously compounded fwd rate at date/time dt."""
    #     t = input_time(dt, self)
    #     epsilon = 1.0e-6
    #     df1 = self.df_t(t)
    #     df2 = self.df_t(t + epsilon)
    #     fwd = np.log(df1 / df2) / epsilon
    #     return fwd

    # ###########################################################################

    # def fwd_rate(self, date1: Date, date2: Date, dc_type: DayCountTypes):
    #     """Calculate the forward rate according to the specified
    #     day count convention."""

    #     if date1 < self.value_dt:
    #         raise FinError("Date1 before curve value date.")

    #     if date2 < date1:
    #         raise FinError("Date2 must not be before Date1")

    #     day_count = DayCount(dc_type)
    #     year_frac = day_count.year_frac(date1, date2)[0]
    #     df1 = self.df(date1)
    #     df2 = self.df(date2)
    #     fwd = (df1 / df2 - 1.0) / year_frac
    #     return fwd

    ###########################################################################

    def plot_zero_rates(self, title: str):
        """Display yield curve."""

        plt.figure(figsize=(12, 6))
        plt.title(title)
        plt.xlabel("Time to Maturity (years)")
        plt.ylabel("Zero Rate (%)")

        tmax = np.max(self.years_to_maturity)
        t = np.linspace(1.0e-6, float(tmax), 100)

        zero_rate = self.zero_rate_t(t)
        zero_rate = scale(zero_rate, 100.0)
        plt.plot(t, zero_rate, label="Zero Rate Bootstrap", marker="o")
        plt.legend(loc="lower right")
        plt.ylim((min(zero_rate) - 0.3, max(zero_rate) * 1.1))
        plt.grid(True)

    ###########################################################################

    def plot_fwd_rates(self, title: str):
        """Display yield curve."""

        plt.figure(figsize=(12, 6))
        plt.title(title)
        plt.xlabel("Time to Maturity (years)")
        plt.ylabel("Forward Rate (%)")

        tmax = np.max(self.years_to_maturity)
        t = np.linspace(1.0e-6, float(tmax), 100)

        fwd_rate = self.fwd_rate_inst_t(t)
        zero_rate = scale(fwd_rate, 100.0)
        plt.plot(t, zero_rate, label="Fwd Rate Bootstrap", marker="o")
        plt.legend(loc="lower right")
        plt.ylim((min(zero_rate) - 0.3, max(zero_rate) * 1.1))
        plt.grid(True)

    ###########################################################################

    def __repr__(self):
        # TODO
        header = "TIMES,DISCOUNT FACTORS"
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        value_table = [self._times, self._dfs]
        precision = "10.7f"
        s += table_to_string(header, value_table, precision)
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
