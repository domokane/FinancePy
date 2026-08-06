# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import os
from typing import Union

import numpy as np

from .interpolator import Interpolator, InterpTypes, interpolate

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.frequency import annual_frequency, FrequencyTypes
from ...utils.compounding import CompoundingTypes
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.math import test_monotonicity
from ...utils.global_vars import G_SMALL
from ...utils.schedule import Schedule
from ...utils.helpers import check_argument_types
from ...utils.helpers import times_from_dates
from ...utils.helpers import label_to_string

import matplotlib.pyplot as plt

########################################################################################


class DiscountCurve:
    """This is a base discount curve which has an internal representation of
    a vector of times and discount factors and an interpolation scheme for
    interpolating between these fixed points."""

    ####################################################################################

    def __init__(
        self,
        value_dt: Date,
        df_dates: list = None,
        df_values: Union[list, np.ndarray] = None,
        interp_type: InterpTypes = InterpTypes.FLAT_FWD_RATES,
        time_dc_type: DayCountTypes = DayCountTypes.ACT_365F,
    ):
        """Create the discount curve from a vector of times and discount
        factors with an anchor date and specify an interpolation scheme. As we
        are explicitly linking dates and discount factors, we do not need to
        specify any compounding convention or day count calculation since
        discount factors are pure prices. We do however need to specify a
        convention for interpolating the discount factors in time."""
        check_argument_types(self.__init__, locals())

        if not isinstance(time_dc_type, DayCountTypes):
            raise FinError("Invalid time day count type.")

        self.time_dc_type = time_dc_type

        # Validate curve
        if df_dates is None:
            df_dates = [value_dt]

        if df_values is None:
            df_values = np.array([1.0], dtype=float)
        else:
            df_values = np.asarray(df_values, dtype=float)

        if len(df_dates) != len(df_values):
            raise FinError("df_dates and df_values must have the same length.")

        if np.any(~np.isfinite(df_values)) or np.any(df_values <= 0.0):
            raise FinError("Discount factors must be finite and positive")

        # The internal representation of times and dfs is hidden but
        # access is controlled using getters and setters
        self._times = [0.0]
        self._dfs = [1.0]
        self._df_dates = [value_dt]

        num_points = len(df_dates)
        start_index = 0

        if num_points > 0 and df_dates[0] == value_dt:
            if np.abs(df_values[0] - 1.0) > 1e-6:
                raise FinError("Value date discount factor should equal 1.0")
            self._dfs[0] = df_values[0]
            start_index = 1

        for i in range(start_index, num_points):
            t = times_from_dates(value_dt, df_dates[i], time_dc_type)
            self._times.append(t)
            self._dfs.append(df_values[i])
            self._df_dates.append(df_dates[i])

        self._times = np.array(self._times)
        self._dfs = np.array(self._dfs)

        if test_monotonicity(self._times) is False:
            print(self._times)
            raise FinError("Times are not sorted in increasing order")

        self.value_dt = value_dt
        self.freq_type = FrequencyTypes.CONTINUOUS

        self._interp_type = interp_type
        self._interpolator = Interpolator(self._interp_type)
        self.fit(self._times, self._dfs)

    ####################################################################################

    def set_times(self, times: np.ndarray):

        self._times = np.asarray(times, dtype=float)

        if not test_monotonicity(self._times):
            raise FinError("Times are not sorted in increasing order")

    ####################################################################################

    def set_dfs(self, dfs: np.ndarray):
        """Set the discount factor at the last maturity time."""

        if np.any(~np.isfinite(dfs)) or np.any(dfs <= 0.0):
            raise FinError("Discount factors must be finite and positive")

        self._dfs = np.asarray(dfs, dtype=float)

    ####################################################################################

    def set_last_df(self, df):
        """Set the discount factor at the last maturity time."""

        if np.any(~np.isfinite(df)) or np.any(df <= 0.0):
            raise FinError("Discount factor must be finite and positive")

        n_points = len(self._dfs)
        self._dfs[n_points - 1] = df

    ####################################################################################

    def fit(self, times: np.ndarray, dfs: np.ndarray):
        """Fit the interpolator to the given times and discount factors."""

        self._interpolator.fit(times, dfs)

    ####################################################################################

    def fwd_rate(
        self,
        start_dt: Union[list, Date],
        date_or_tenor: Union[Date, str, list],
        accrual_dc_type: DayCountTypes = DayCountTypes.ACT_360,
        comp_type: CompoundingTypes = CompoundingTypes.MMKT,
        freq_type: FrequencyTypes = None,
    ):
        """Simple-compounded forward rate between dates.

        accrual_dc_type is used for the forward accrual denominator.
        self.time_dc_type is used for curve time conversion / interpolation.
        """

        if not isinstance(accrual_dc_type, DayCountTypes):
            raise FinError("Invalid Accrual Day Count type.")

        if not isinstance(comp_type, CompoundingTypes):
            raise FinError("Invalid Compounding type.")

        if comp_type == CompoundingTypes.COMPOUNDED and freq_type == None:
            raise FinError("Require a Frequency type for COMPOUNDING.")

        scalar_input = isinstance(start_dt, Date)
        start_dts = [start_dt] if scalar_input else start_dt

        if not isinstance(start_dts, list):
            raise FinError("Start date must be a Date or list of Dates.")

        if len(start_dts) == 0:
            return np.array([])

        if not all(isinstance(dt, Date) for dt in start_dts):
            raise FinError("start_dt list must contain only Dates.")

        if isinstance(date_or_tenor, str):
            end_dts = [dt.add_tenor(date_or_tenor) for dt in start_dts]

        elif isinstance(date_or_tenor, Date):
            end_dts = [date_or_tenor for _ in start_dts]

        elif isinstance(date_or_tenor, list):
            if len(date_or_tenor) != len(start_dts):
                raise FinError("date_or_tenor list must match start_dt list length.")

            if not all(isinstance(dt, Date) for dt in date_or_tenor):
                raise FinError("date_or_tenor list must contain only Dates.")

            end_dts = date_or_tenor

        else:
            raise FinError(
                "date_or_tenor must be a Date, tenor string, or list of Dates."
            )

        day_counter = DayCount(accrual_dc_type)

        accruals = np.array(
            [
                day_counter.year_frac(dt1, dt2)[0]
                for dt1, dt2 in zip(start_dts, end_dts)
            ],
            dtype=float,
        )

        if np.any(accruals <= 0.0):
            raise FinError("Forward end date must be after start date.")

        t1 = times_from_dates(self.value_dt, start_dts, self.time_dc_type)
        t2 = times_from_dates(self.value_dt, end_dts, self.time_dc_type)

        rates = self.fwd_rate_t(t1, t2, accruals, comp_type, freq_type)

        if scalar_input:
            return float(rates[0])
        else:
            return rates

    ###########################################################################

    def fwd_rate_t(
        self,
        start_t: Union[float, list, np.ndarray],
        end_t: Union[float, list, np.ndarray],
        accrual: Union[float, list, np.ndarray] = None,
        comp_type: CompoundingTypes = CompoundingTypes.MMKT,
        freq_type: FrequencyTypes = None,
    ):
        """Forward rate between curve times."""

        if comp_type == CompoundingTypes.COMPOUNDED and freq_type is None:
            raise FinError("Require Frequency type for COMPOUNDED interest.")

        t1, scalar_input_1 = self._to_time_array(start_t)
        t2, scalar_input_2 = self._to_time_array(end_t)

        if len(t1) != len(t2):
            if len(t1) == 1:
                t1 = np.full_like(t2, t1[0])
            elif len(t2) == 1:
                t2 = np.full_like(t1, t2[0])
            else:
                raise FinError("Start and end time vectors have different lengths.")

        if accrual is None:
            accruals = t2 - t1
        else:
            accruals = np.array(accrual, ndmin=1, dtype=float)

            if len(accruals) != len(t1):
                if len(accruals) == 1:
                    accruals = np.full_like(t1, accruals[0])
                else:
                    raise FinError("Accrual vector has incorrect length.")

        if np.any(t2 <= t1):
            raise FinError("End time must be after start time.")

        if np.any(accruals <= 0.0):
            raise FinError("Accrual factor must be positive.")

        df1 = self.df_t(t1)
        df2 = self.df_t(t2)
        df_ratio = df1 / df2

        if comp_type == CompoundingTypes.CONTINUOUS:
            fwd_rates = np.log(df_ratio) / accruals

        elif comp_type == CompoundingTypes.MMKT:
            fwd_rates = (df_ratio - 1.0) / accruals

        else:
            f = annual_frequency(freq_type)
            fwd_rates = (np.power(df_ratio, 1.0 / (f * accruals)) - 1.0) * f

        scalar_output = scalar_input_1 and scalar_input_2
        if scalar_output:
            return float(fwd_rates[0])
        else:
            return fwd_rates

    ###########################################################################

    def fwd_rate_inst(self, dts: Union[Date, list], dt: float = 1.0e-6):
        """Instantaneous continuously compounded forward rate at date(s)."""
        times = times_from_dates(self.value_dt, dts, self.time_dc_type)
        return self.fwd_rate_inst_t(times, dt)

    ###########################################################################

    def fwd_rate_inst_t(
        self,
        t: Union[float, list, np.ndarray],
        dt: float = 1.0e-4,
    ):
        """Instantaneous continuously compounded forward rate at time t.

        Computes:

            f(t) = - d ln Z(t) / dt

        using a central finite difference.
        """

        times, scalar_input = self._to_time_array(t)
        times = np.maximum(times, dt)

        df_minus = self.df_t(times - dt)
        df_plus = self.df_t(times + dt)

        fwds = np.log(df_minus / df_plus) / (2.0 * dt)

        if scalar_input:
            return float(fwds[0])
        else:
            return fwds

    ###########################################################################

    def zero_rate(
        self,
        maturity_dt: Union[list, Date],
        freq_type: FrequencyTypes = FrequencyTypes.ANNUAL,
    ):
        """Calculate zero rates with specified compounding frequency."""

        if not isinstance(freq_type, FrequencyTypes):
            raise FinError("Invalid Frequency type.")

        times = times_from_dates(self.value_dt, maturity_dt, self.time_dc_type)
        zero_rates = self.zero_rate_t(times, freq_type)
        return zero_rates

    ###########################################################################

    def zero_rate_cc(self, maturity_dt: Union[list, Date]):
        """Calculate zero rates with continuous compounding."""

        times = times_from_dates(self.value_dt, maturity_dt, self.time_dc_type)
        zero_rates_cc = self.zero_rate_t(times, FrequencyTypes.CONTINUOUS)
        return zero_rates_cc

    ####################################################################################

    def zero_rate_cc_t(self, t):
        """Calculate zero rates with continuous compounding."""
        zero_rates_cc = self.zero_rate_t(t, FrequencyTypes.CONTINUOUS)
        return zero_rates_cc

    ####################################################################################

    def zero_rate_t(
        self,
        t,
        freq_type: FrequencyTypes = FrequencyTypes.ANNUAL,
    ):
        if not isinstance(freq_type, FrequencyTypes):
            raise FinError("Invalid Frequency type.")

        if np.any(t < 0):
            raise FinError("Times to maturity cannot be negative.")

        times, scalar_input = self._to_time_array(t)
        times = np.maximum(times, G_SMALL)

        dfs = self.df_t(times)
        f = annual_frequency(freq_type)

        if freq_type == FrequencyTypes.CONTINUOUS:
            rates = -np.log(dfs) / times
        elif freq_type == FrequencyTypes.SIMPLE_INTEREST:
            rates = (1.0 / dfs - 1.0) / times
        else:
            rates = (np.power(dfs, -1.0 / (times * f)) - 1.0) * f

        return float(rates[0]) if scalar_input else rates

    ####################################################################################

    def swap_rate(
        self,
        effective_dt: Date,
        maturity_dt: Union[list, Date],
        freq_type: FrequencyTypes = FrequencyTypes.ANNUAL,
        accrual_dc_type: DayCountTypes = DayCountTypes.THIRTY_E_360,
    ):
        rate = self.par_rate(effective_dt, maturity_dt, freq_type, accrual_dc_type)
        return rate

    ####################################################################################

    def swap_rate_old(
        self,
        effective_dt: Date,
        maturity_dt: Union[list, Date],
        freq_type: FrequencyTypes = FrequencyTypes.ANNUAL,
        accrual_dc_type: DayCountTypes = DayCountTypes.THIRTY_E_360,
    ):
        """Calculate the swap rate to maturity date. This is the rate paid by
        a swap that has a price of par today. This is the same as a Libor swap
        rate except that we do not do any business day adjustments."""

        # Note that this function does not call the IborSwap class to
        # calculate the swap rate since that will create a circular dependency.
        # I therefore recreate the actual calculation of the swap rate here.

        if effective_dt < self.value_dt:
            raise FinError("Swap starts before the curve valuation date.")

        if isinstance(freq_type, FrequencyTypes) is False:
            raise FinError("Invalid Frequency type.")

        if not isinstance(accrual_dc_type, DayCountTypes):
            raise FinError("Invalid Day Count type.")

        if freq_type == FrequencyTypes.SIMPLE_INTEREST:
            raise FinError("Cannot calculate par rate with simple yield freq.")

        if freq_type == FrequencyTypes.CONTINUOUS:
            raise FinError("Cannot calculate par rate with continuous freq.")

        scalar_input = isinstance(maturity_dt, Date)

        if scalar_input:
            maturity_dts = [maturity_dt]
        else:
            maturity_dts = maturity_dt

        df_start = self.df(effective_dt)
        accrual_day_counter = DayCount(accrual_dc_type)
        par_rates = []

        for maturity_dt in maturity_dts:

            schedule = Schedule(effective_dt, maturity_dt, freq_type)

            flow_dts = schedule.generate()
            flow_dts[0] = effective_dt

            pv01 = 0.0
            df_end = None
            prev_dt = effective_dt

            for next_dt in flow_dts[1:]:
                df_end = self.df(next_dt)
                alpha = accrual_day_counter.year_frac(prev_dt, next_dt)[0]
                pv01 += alpha * df_end
                prev_dt = next_dt

            if abs(pv01) < G_SMALL:
                par_rate = 0.0
            else:
                par_rate = (df_start - df_end) / pv01

            par_rates.append(par_rate)

        par_rates = np.array(par_rates)

        if scalar_input:
            return par_rates[0]
        else:
            return par_rates

    ###########################################################################

    def par_rate(
        self,
        effective_dt: Date,
        maturity_dt: Union[list, Date],
        freq_type: FrequencyTypes = FrequencyTypes.ANNUAL,
        accrual_dc_type: DayCountTypes = DayCountTypes.THIRTY_E_360,
    ):
        """Calculate the swap rate to maturity date. This is the rate paid by
        a swap that has a price of par today. This is the same as a Libor swap
        rate except that we do not do any business day adjustments."""

        # Note that this function does not call the IborSwap class to
        # calculate the swap rate since that will create a circular dependency.
        # I therefore recreate the actual calculation of the swap rate here.

        if effective_dt < self.value_dt:
            raise FinError("Swap starts before the curve valuation date.")

        if isinstance(freq_type, FrequencyTypes) is False:
            raise FinError("Invalid Frequency type.")

        if not isinstance(accrual_dc_type, DayCountTypes):
            raise FinError("Invalid Day Count type.")

        if freq_type == FrequencyTypes.SIMPLE_INTEREST:
            raise FinError("Cannot calculate par rate with simple yield freq.")

        if freq_type == FrequencyTypes.CONTINUOUS:
            raise FinError("Cannot calculate par rate with continuous freq.")

        scalar_input = isinstance(maturity_dt, Date)

        if scalar_input:
            maturity_dts = [maturity_dt]
        else:
            maturity_dts = maturity_dt

        acc_day_counter = DayCount(accrual_dc_type)

        t_start = times_from_dates(self.value_dt, effective_dt, self.time_dc_type)

        par_rates = []

        for mat_dt in maturity_dts:

            if mat_dt <= effective_dt:
                raise FinError("Swap maturity date must be after effective date.")

            schedule = Schedule(effective_dt, mat_dt, freq_type)
            flow_dts = schedule.generate()

            payment_times = np.array(
                [
                    times_from_dates(self.value_dt, dt, self.time_dc_type)
                    for dt in flow_dts[1:]
                ]
            )

            accrual_factors = []
            for prev_dt, next_dt in zip(flow_dts[:-1], flow_dts[1:]):
                accd = acc_day_counter.year_frac(prev_dt, next_dt)[0]
                accrual_factors.append(accd)

            accrual_factors = np.array(accrual_factors)

            rate = self.par_rate_t(t_start, payment_times, accrual_factors)
            par_rates.append(rate)

        par_rates = np.array(par_rates)

        if scalar_input:
            return par_rates[0]

        return par_rates

    ###########################################################################

    def par_rate_simple_t(self, t_maturity, freq=2):
        """
        Simple spot-starting par rate from curve times. Handles bonds with coupons
        equal to c/f. Calculates payments backwards from maturity. Not exact for swaps
        as it ignores accrual factor weighted coupon payments.

        Parameters
        ----------
        maturity : float
            Final maturity in years.
        freq : int
            Fixed-leg payments per year.
        """

        dt = 1.0 / freq
        payment_times = [t_maturity]
        t = t_maturity

        while t > dt:
            t -= dt
            payment_times.append(t)

        payment_times = np.array(sorted(payment_times), dtype=float)
        accrual_factors = np.diff(np.concatenate(([0.0], payment_times)))
        p = self.par_rate_t(0.0, payment_times, accrual_factors)
        return p

    #############################################################################

    def par_rate_t(
        self,
        t_start: float,
        payment_times: np.ndarray,
        accrual_factors: np.ndarray,
    ):

        df_start = self.df_t(t_start)
        dfs = self.df_t(payment_times)

        pv01 = np.sum(accrual_factors * dfs)

        if abs(pv01) < G_SMALL:
            return 0.0

        swap_rate = (df_start - dfs[-1]) / pv01
        return swap_rate

    ###########################################################################

    def df(self, dt: Union[list, Date]):
        """Function to calculate a discount factor from a date or a
        vector of dates. The time day count determines how dates get converted
        to years."""

        times = times_from_dates(self.value_dt, dt, self.time_dc_type)
        dfs = self.df_t(times)

        if isinstance(dt, Date):
            if np.isscalar(dfs):
                return float(dfs)
            else:
                return float(dfs[0])

        return np.asarray(dfs, dtype=float)

    ###########################################################################

    def df_t(self, t: Union[float, np.ndarray]):
        """Function to calculate a discount factor from a time or a
        vector of times. Discourage usage in favour of passing in dates."""

        times, scalar_input = self._to_time_array(t)

        if self._interp_type in (
            InterpTypes.FLAT_FWD_RATES,
            InterpTypes.LINEAR_DISCOUNT,
            InterpTypes.LINEAR_ZERO_RATES,
        ):
            dfs = interpolate(times, self._times, self._dfs, self._interp_type.value)
        else:
            dfs = self._interpolator.interpolate(times)

        if scalar_input:
            return float(dfs[0])
        else:
            return dfs

    ###########################################################################

    def bump_parallel(self, bump_size: float):
        """Adjust the continuously compounded forward rates by a perturbation
        upward equal to the bump size and return a curve object with this bumped
        curve. This is used for interest rate risk."""

        bumped_dfs = self._dfs * np.exp(-bump_size * self._times)
        bumped_dfs[0] = self._dfs[0]

        bumped_disc_curve = DiscountCurve(
            self.value_dt,
            self._df_dates,
            bumped_dfs,
            self._interp_type,
            self.time_dc_type,
        )

        return bumped_disc_curve

    ###########################################################################

    def bump_bucket(
        self,
        bucket_start: float,  # start in years
        bucket_end: float,  # end in years
        bump_size: float,  # e.g. 0.0001 = +1bp
        interp_type=None,
    ):
        """
        Apply a bucket (key rate) shift: shift continuous forward rates
        only between bucket_start and bucket_end years.
        """

        if interp_type is None:
            interp_type = self._interp_type

        if bucket_start < 0.0:
            raise FinError("Bucket start must be non-negative.")

        if bucket_end <= bucket_start:
            raise FinError("Bucket end must be after bucket start.")

        times = self._times.copy()
        dfs = self._dfs.copy()

        shifted_lengths = np.maximum(0.0, np.minimum(times, bucket_end) - bucket_start)

        dfs = dfs * np.exp(-bump_size * shifted_lengths)

        # Build new curve
        bumped_curve = DiscountCurve(
            value_dt=self.value_dt,
            df_dates=self._df_dates,
            df_values=dfs,
            interp_type=interp_type,
            time_dc_type=self.time_dc_type,
        )

        return bumped_curve

    ###########################################################################

    def _zero_to_df(
        self,
        rates: Union[float, np.ndarray],
        times: Union[float, np.ndarray],
        freq_type: FrequencyTypes,
    ):
        """Convert a zero with a specified compounding frequency and day count
        convention to a discount factor for a single maturity date or a list of
        dates."""

        times, scalar_input = self._to_time_array(times)
        rates = np.asarray(rates, dtype=float)

        t = np.maximum(times, G_SMALL)
        f = annual_frequency(freq_type)

        if freq_type == FrequencyTypes.CONTINUOUS:
            dfs = np.exp(-rates * t)
        elif freq_type == FrequencyTypes.SIMPLE_INTEREST:
            dfs = 1.0 / (1.0 + rates * t)
        elif freq_type in {
            FrequencyTypes.ANNUAL,
            FrequencyTypes.SEMI_ANNUAL,
            FrequencyTypes.QUARTERLY,
            FrequencyTypes.MONTHLY,
        }:
            dfs = 1.0 / np.power(1.0 + rates / f, f * t)
        else:
            raise FinError("Unknown Frequency type")

        return np.asarray(dfs, dtype=float)

    ###########################################################################

    def survival_prob(self, dt: Date) -> float:
        """This returns a survival probability to a specified date based on
        the assumption that the continuously compounded rate is a default
        hazard rate in which case the survival probability is directly
        analogous to a discount factor."""

        q = self.df(dt)
        return q

    ###########################################################################

    def _to_time_array(self, t):
        """Convert a scalar/list/array time input into a 1D NumPy array."""

        scalar_input = np.isscalar(t)

        if scalar_input:
            return np.array([float(t)], dtype=float), True

        return np.asarray(t, dtype=float), False

    ###########################################################################

    def plot(
        self,
        title,
        times: np.ndarray = None,
        ymin: float = None,
        ymax: float = None,
        filename: str = None,
    ):
        """Display discount curve."""

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
        plt.title(title)

        if times is None:
            tmax = np.max(self._times)
            times = np.linspace(0.0, tmax, int(tmax * 12.0))
        else:
            times = np.asarray(times, dtype=float)

        if np.any(times < 0.0):
            raise FinError("Plot times must be strictly positive.")

        times = np.maximum(times, 1e-8)

        zeros = self.zero_rate_cc_t(times)
        fwds = self.fwd_rate_inst_t(times)

        # Bump the forwards by 1bp in case the curve is flat so they can be seen
        fwds = fwds + 1e-4

        plt.plot(times, zeros * 100, label="CC Zero Rates")
        plt.plot(times, fwds * 100, label="Inst Fwd Rates")

        plt.xlabel("Time to Maturity (years)")
        plt.ylabel("Rate (%)")
        plt.legend()

        if ymin is not None and ymax is not None:
            plt.ylim(ymin, ymax)

        plt.grid(True, alpha=0.3)
        plt.tight_layout()

        if filename is not None:
            plt.savefig(filename, bbox_inches="tight", pad_inches=0.02)

        plt.show()

    ###########################################################################

    def __repr__(self):

        # Hardcode this as we want this not parent class
        s = label_to_string("OBJECT_TYPE", "DiscountCurve")
        s += label_to_string("VALUE DATE", (self.value_dt))

        s += "    DATES      TIMES(YRS) DISC FACTORS\n"
        for dt, t, df in zip(self._df_dates, self._times, self._dfs):
            s += label_to_string(f"{str(dt):>12}", f"{t:12.6f}", f"{df:12.8f}\n")

        if hasattr(self, "_interp_type") and self._interp_type is not None:
            s += label_to_string("INTERPOLATION TYPE", self._interp_type.name)
        else:
            s += label_to_string("INTERPOLATION TYPE", "N/A")

        s += label_to_string("TIME DAY COUNT TYPE", (self.time_dc_type.name))

        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)
