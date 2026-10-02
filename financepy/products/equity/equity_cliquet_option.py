##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...utils.frequency import FrequencyTypes
from ...utils.error import FinError
from ...utils.global_types import OptionTypes
from ...utils.global_types import CliquetTypes
from ...utils.global_types import PaymentTimingTypes

from ...utils.helpers import label_to_string, check_argument_types
from ...utils.helpers import option_years
from ...utils.date import Date
from ...utils.calendar import BusDayAdjustTypes
from ...utils.calendar import CalendarTypes, DateGenRuleTypes
from ...utils.schedule import Schedule
from ...utils.check_values import check_curve_dt
from ...utils.stats import std_err
from ...utils.mc_result import MCResult

from ...products.equity.equity_option import EquityOption
from ...market.curves.flat_discount_curve import DiscountCurve

from ...models.black_scholes_analytic import european_value
from ...models.black_scholes import BlackScholes
from ...models.model import Model

########################################################################################
# TODO: Do we need to day count adjust option payoffs ?
# TODO: Monte Carlo pricer
########################################################################################


class EquityCliquetOption(EquityOption):
    """A EquityCliquetOption is a series of options which start and stop at
    successive times with each subsequent option resetting its strike to be ATM
    at the start of its life. This is also known as a reset option."""

    def __init__(
        self,
        start_dt: Date,
        final_expiry_dt: Date,
        opt_type: OptionTypes,
        freq_type: FrequencyTypes,
        payoff_type: CliquetTypes = CliquetTypes.PRICE,
        payoff_timing_type: PaymentTimingTypes = PaymentTimingTypes.PERIODIC,
        notional: float = 1.0,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
    ) -> None:
        """Create the EquityCliquetOption by passing in the start date
        and the end date and whether it is a call or a put. Some additional
        data is needed in order to calculate the individual payments."""

        check_argument_types(self.__init__, locals())

        if opt_type != OptionTypes.EUROPEAN_CALL and opt_type != OptionTypes.EUROPEAN_PUT:
            raise FinError("Unknown OPTION_TYPE" + str(opt_type))

        if final_expiry_dt < start_dt:
            raise FinError("Expiry date precedes start date")

        self.start_dt = start_dt
        self.final_expiry_dt = final_expiry_dt
        self.opt_type = opt_type
        self.freq_type = freq_type
        self.payoff_type = payoff_type
        self.payoff_timing_type = payoff_timing_type
        self.notional = notional
        self.cal_type = cal_type
        self.bd_type = bd_type
        self.dg_type = dg_type

        self._validate()

        self.v_options = None
        self._dfs = None
        self.actual_dts = None

        self.expiry_dts = Schedule(
            self.start_dt,
            self.final_expiry_dt,
            self.freq_type,
            self.cal_type,
            self.bd_type,
            self.dg_type,
        ).generate()

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
    ):
        """Value the cliquet option as a sequence of options using the Black-
        Scholes model."""

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.final_expiry_dt:
            raise FinError("Value date after final expiry date.")

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        s0 = stock_price
        v_cliquet = 0.0

        self.v_options = []
        self.dfs = []
        self.actual_dts = []

        if isinstance(model, BlackScholes):

            fwd_vol = model.volatility
            fwd_vol = max(fwd_vol, 1e-6)

        else:
            raise FinError("Unknown Model Type")

        df_final = discount_curve.df(self.final_expiry_dt)

        dt_prev = value_dt

        for dt in self.expiry_dts:

            if dt > value_dt:

                t_vol = option_years(dt_prev, dt)

                df_start = discount_curve.df(dt_prev)
                df_end = discount_curve.df(dt)
                dq_start = dividend_curve.df(dt_prev)

                fwd_r = discount_curve.fwd_zero_rate_cc(dt_prev, dt)
                fwd_q = dividend_curve.fwd_zero_rate_cc(dt_prev, dt)

                v = european_value(
                    1.0,
                    t_vol,
                    1.0,
                    fwd_r,
                    fwd_q,
                    fwd_vol,
                    self.opt_type.value,
                )

                # Value today assuming payment at period expiry dt
                if self.payoff_type == CliquetTypes.PRICE:
                    period_value = s0 * dq_start * v

                elif self.payoff_type == CliquetTypes.RETURN:
                    period_value = df_start * v

                else:
                    raise FinError(
                        "Unknown CLIQUET_PAYOFF_TYPE "
                        + str(self.payoff_type)
                    )

                # If the payoff is deferred from dt to final maturity,
                # discount the locked amount over [dt, final_expiry_dt].
                if self.payoff_timing_type == PaymentTimingTypes.MATURITY:

                    period_value *= df_final / df_end

                elif self.payoff_timing_type == PaymentTimingTypes.PERIODIC:

                    pass

                else:

                    raise FinError(
                        "Unknown PAYMENT_TIMING_TYPE "
                        + str(self.payoff_timing_type)
                    )

                v_cliquet += self.notional * period_value

                self.dfs.append(df_end)
                self.v_options.append(period_value)
                self.actual_dts.append(dt)

                dt_prev = dt

        return v_cliquet

    ##########################################################################

    def value_mc(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        num_paths: int = 200_000,
        seed: int = 4242,
    ):
        """Value the cliquet option using Monte Carlo simulation.

        The stock is simulated exactly under Black-Scholes from reset date
        to reset date, so there is no Euler time-discretisation error.

        Returns
        -------
        value : float
            Monte Carlo value.

        stderr : float
            Monte Carlo standard error estimated from antithetic pair
            averages.

        Notes
        -----
        This implementation requires value_dt == start_dt. Valuation after
        inception requires historical reset fixings and, for MATURITY
        payment timing, already-realised period payoffs.
        """

        if not isinstance(value_dt, Date):
            raise FinError("Valuation date is not a Date")

        if value_dt != self.start_dt:
            raise FinError(
                "Monte Carlo valuation requires value_dt == start_dt."
            )

        if value_dt > self.final_expiry_dt:
            raise FinError("Value date after final expiry date.")

        if stock_price <= 0.0:
            raise FinError("Stock price must be positive")

        if num_paths < 2:
            raise FinError("Number of paths must be at least 2")

        if not isinstance(model, BlackScholes):
            raise FinError(
                "Monte Carlo valuation requires a BlackScholes model"
            )

        if self.opt_type != OptionTypes.EUROPEAN_CALL:
            raise FinError(
                "Monte Carlo valuation currently supports call cliquets only"
            )

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        # Force an even number of paths so that every path has
        # an antithetic partner.
        num_half_paths = (num_paths + 1) // 2
        num_paths = 2 * num_half_paths

        sigma = max(model.volatility, 1.0e-12)

        rng = np.random.default_rng(seed)

        # Keep antithetic paths separate so that the standard error can
        # be calculated from independent pair averages.
        s_plus = np.full(
            num_half_paths,
            stock_price,
            dtype=float,
        )

        s_minus = np.full(
            num_half_paths,
            stock_price,
            dtype=float,
        )

        pv_plus = np.zeros(num_half_paths)
        pv_minus = np.zeros(num_half_paths)

        df_final = discount_curve.df(
            self.final_expiry_dt
        )

        dt_prev = value_dt

        # -------------------------------------------------------------
        # Simulate each reset period
        # -------------------------------------------------------------

        for dt in self.expiry_dts:

            if dt <= value_dt:
                continue

            t = option_years(
                dt_prev,
                dt,
            )

            if t <= 0.0:
                dt_prev = dt
                continue

            # Forward rates applying over this reset period.
            fwd_r = discount_curve.fwd_zero_rate_cc(
                dt_prev,
                dt,
            )

            fwd_q = dividend_curve.fwd_zero_rate_cc(
                dt_prev,
                dt,
            )

            # Independent normals for the antithetic pairs.
            z = rng.standard_normal(
                num_half_paths
            )

            drift = (
                fwd_r
                - fwd_q
                - 0.5 * sigma * sigma
            ) * t

            diffusion = (
                sigma
                * np.sqrt(t)
                * z
            )

            # Exact GBM transitions under the risk-neutral measure.
            s_next_plus = (
                s_plus
                * np.exp(drift + diffusion)
            )

            s_next_minus = (
                s_minus
                * np.exp(drift - diffusion)
            )

            # ---------------------------------------------------------
            # Period payoff
            # ---------------------------------------------------------

            if self.payoff_type == CliquetTypes.PRICE:

                payoff_plus = np.maximum(
                    s_next_plus - s_plus,
                    0.0,
                )

                payoff_minus = np.maximum(
                    s_next_minus - s_minus,
                    0.0,
                )

            elif self.payoff_type == CliquetTypes.RETURN:

                ret_plus = (
                    s_next_plus / s_plus - 1.0
                )

                ret_minus = (
                    s_next_minus / s_minus - 1.0
                )

                payoff_plus = np.maximum(
                    ret_plus,
                    0.0,
                )

                payoff_minus = np.maximum(
                    ret_minus,
                    0.0,
                )

            else:

                raise FinError(
                    "Unknown CLIQUET_PAYOFF_TYPE "
                    + str(self.payoff_type)
                )

            # ---------------------------------------------------------
            # Payment timing
            # ---------------------------------------------------------

            if (
                self.payoff_timing_type
                == PaymentTimingTypes.PERIODIC
            ):

                df_pay = discount_curve.df(dt)

            elif (
                self.payoff_timing_type
                == PaymentTimingTypes.MATURITY
            ):

                df_pay = df_final

            else:

                raise FinError(
                    "Unknown PAYMENT_TIMING_TYPE "
                    + str(self.payoff_timing_type)
                )

            pv_plus += (
                df_pay
                * payoff_plus
            )

            pv_minus += (
                df_pay
                * payoff_minus
            )

            # The ending stock becomes the reset strike for
            # the next period.
            s_plus = s_next_plus
            s_minus = s_next_minus

            dt_prev = dt

        # -------------------------------------------------------------
        # Antithetic estimator
        # -------------------------------------------------------------

        pair_values = 0.5 * (
            pv_plus + pv_minus
        )

        pair_values *= self.notional

        value = np.mean(pair_values)

        e = std_err(pair_values)

        return MCResult(value, e)

    ##########################################################################

    def _validate(self):
        """Validate EquityCliquetOption constructor inputs."""

        # Dates
        if not isinstance(self.start_dt, Date):
            raise FinError("start_dt must be a Date")

        if not isinstance(self.final_expiry_dt, Date):
            raise FinError("final_expiry_dt must be a Date")

        if self.final_expiry_dt <= self.start_dt:
            raise FinError("final_expiry_dt must be after start_dt")

        # Option type
        if self.opt_type not in (
            OptionTypes.EUROPEAN_CALL,
            OptionTypes.EUROPEAN_PUT,
        ):
            raise FinError(
                "opt_type must be EUROPEAN_CALL or EUROPEAN_PUT"
            )

    ###########################################################################

    def print_payments(self):

        if self.v_options is None:
            raise FinError("Options not created yet.")

        num_options = len(self.v_options)
        for i in range(0, num_options):
            print(self.actual_dts[i], self.dfs[i], self.v_options[i])

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("START_DATE", self.start_dt)
        s += label_to_string("FINAL EXPIRY DATE", self.final_expiry_dt)
        s += label_to_string("OPTION_TYPE", self.opt_type)
        s += label_to_string("FREQUENCY_TYPE", self.freq_type)
        s += label_to_string("PAYOFF_TYPE", self.payoff_type)
        s += label_to_string("PAYMENT_TIME_TYPE", self.payoff_timing_type)
        s += label_to_string("CALENDAR_TYPE", self.cal_type)
        s += label_to_string("BUS_DAY_ADJUST_TYPE", self.bd_type)
        s += label_to_string("DATE_GEN_RULE_TYPE", self.dg_type, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
