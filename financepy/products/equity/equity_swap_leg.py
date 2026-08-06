##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from typing import Union

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.calendar import CalendarTypes, DateGenRuleTypes
from ...utils.calendar import Calendar, BusDayAdjustTypes
from ...utils.schedule import Schedule
from ...utils.helpers import (
    format_table,
    label_to_string,
    check_argument_types,
)
from ...utils.global_types import SwapTypes, ReturnTypes
from ...market.curves.discount_curve import DiscountCurve
from ...market.curves.flat_discount_curve import FlatDiscountCurve

##########################################################################


class EquitySwapLeg:
    """Class for managing the equity leg of an equity swap. An equity leg is
    a leg with a sequence of flows calculated according to an ISDA schedule
    and follows the economics of a collection of equity forward contracts.
    """

    def __init__(
        self,
        effective_dt: Date,  # Contract starts or last equity reset
        term_dt_or_tenor: Union[Date, str],  # Date contract ends
        leg_type: SwapTypes,
        freq_type: FrequencyTypes,
        accrual_dc_type: DayCountTypes,
        strike: float,  # Price at effective date
        quantity: float = 1.0,  # Quantity at effective date
        payment_lag: int = 0,
        return_type: ReturnTypes = ReturnTypes.TOTAL_RETURN,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
        end_of_month: bool = False,
    ):
        """Create the equity leg of a swap contract giving the contract start
        date, its maturity, underlying strike price and quantity, payment
        frequency, day count convention, return type, and other details"""

        check_argument_types(self.__init__, locals())

        if isinstance(term_dt_or_tenor, Date):
            termination_dt = term_dt_or_tenor
        else:
            termination_dt = effective_dt.add_tenor(term_dt_or_tenor)

        calendar = Calendar(cal_type)
        self.maturity_dt = calendar.adjust(termination_dt, bd_type)

        if effective_dt > self.maturity_dt:
            raise FinError("Effective date after maturity date")

        if quantity < 0:
            # Long/Short is defined by leg_type
            raise FinError("Quantity must be non-negative")

        if return_type != ReturnTypes.TOTAL_RETURN:
            raise NotImplementedError("Return Type still not implemented")

        # To generate ISDA pmnt schedule properly we can't use these types
        if freq_type in (
            FrequencyTypes.CONTINUOUS,
            FrequencyTypes.SIMPLE_INTEREST,
        ):
            print(freq_type)
            raise FinError(
                "Cannot generate payment schedule for this frequency type"
            )

        self.effective_dt = effective_dt
        self.termination_dt = termination_dt
        self.leg_type = leg_type
        self.freq_type = freq_type
        self.payment_lag = payment_lag
        self.strike = strike
        self.quantity = quantity
        self.notional = strike * quantity
        self.return_type = return_type

        self.accrual_dc_type = accrual_dc_type
        self.cal_type = cal_type
        self.bd_type = bd_type
        self.dg_type = dg_type
        self.end_of_month = end_of_month

        self.start_accd_dts = []
        self.end_accd_dts = []
        self.payment_dts = []
        self.payment_amounts = []
        self.year_fracs = []
        self.accrued_days = []
        self.rates = []

        self.fwd_rates = []
        self.div_fwd_rates = []
        self.eq_fwd_rates = []
        self.last_notionals = []
        self.payment_dfs = []
        self.payment_pvs = []
        self.cumulative_pvs = []
        self.current_price = None

        self.generate_payment_dts()

    ###########################################################################

    def generate_payment_dts(self):
        """Generate the Equity leg payment dates and accrual factors. Similar
        to swap float leg, payment values can't be generated, as we do not have
        index curve, dividend curve and equity price."""

        schedule = Schedule(
            self.effective_dt,
            self.maturity_dt,
            self.freq_type,
            self.cal_type,
            self.bd_type,
            self.dg_type,
            end_of_month=self.end_of_month,
        )

        schedule_dts = schedule.adjusted_dts

        if len(schedule_dts) < 2:
            raise FinError("Schedule has none or only one date")

        self.start_accd_dts = []
        self.end_accd_dts = []
        self.payment_dts = []
        self.year_fracs = []
        self.accrued_days = []

        prev_dt = schedule_dts[0]

        day_counter = DayCount(self.accrual_dc_type)
        calendar = Calendar(self.cal_type)

        # All the lists end up with the same length
        for next_dt in schedule_dts[1:]:

            self.start_accd_dts.append(prev_dt)
            self.end_accd_dts.append(next_dt)

            if self.payment_lag == 0:
                payment_dt = next_dt
            else:
                payment_dt = calendar.add_business_days(
                    next_dt, self.payment_lag
                )

            self.payment_dts.append(payment_dt)

            (year_frac, num, _) = day_counter.year_frac(prev_dt, next_dt)

            self.year_fracs.append(year_frac)
            self.accrued_days.append(num)

            prev_dt = next_dt

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        discount_curve: DiscountCurve,
        index_curve: DiscountCurve,
        dividend_curve: DiscountCurve = None,
        current_price: float = None
    ):
        """Value the equity leg with payments from an equity price, quantity,
        an index curve and an [optional] dividend curve. Discounting is based
        on a supplied discount curve as of the valuation date supplied.

        Each reset period is treated as an independent forward contract:

            PV_i = Q * S_0 * df_idx(0, t_{i-1}) * df_div(0, t_{i-1}) * df(0, t_i)
                    * [1 / (df_idx(0, t_i) * df_div(0, t_i)) - 1]

        For the current (in-progress) period where t_{i-1} < value_dt <= t_i,
        the period-start equity price is known (current_price or strike), so
        we value it as a single forward from value_dt to t_i.
        """

        if discount_curve is None:
            raise FinError("Discount curve not provided!")

        if index_curve is None:
            index_curve = discount_curve

        if dividend_curve is None:
            dividend_curve = FlatDiscountCurve(value_dt, 0.0)

        if discount_curve.value_dt != value_dt:
            raise FinError("Discount Curve valuation date not same as value date")

        self.current_price = current_price if current_price is not None else self.strike

        # Reset all result lists
        self.fwd_rates = []
        self.div_fwd_rates = []
        self.eq_fwd_rates = []
        self.last_notionals = []
        self.payment_amounts = []
        self.payment_dfs = []
        self.payment_pvs = []
        self.cumulative_pvs = []

        df_value = discount_curve.df(value_dt)
        leg_pv = 0.0
        num_payments = len(self.payment_dts)

        for i_pmnt in range(num_payments):

            payment_dt = self.payment_dts[i_pmnt]
            start_accrued_dt = self.start_accd_dts[i_pmnt]
            end_accrued_dt = self.end_accd_dts[i_pmnt]

            if payment_dt <= value_dt:
                # Period already paid — no value
                self.fwd_rates.append(0.0)
                self.div_fwd_rates.append(0.0)
                self.eq_fwd_rates.append(0.0)
                self.last_notionals.append(0.0)
                self.payment_amounts.append(0.0)
                self.payment_dfs.append(0.0)
                self.payment_pvs.append(0.0)
                self.cumulative_pvs.append(leg_pv)
                continue

            # ------------------------------------------------------------------
            # Discount factors from the curve's anchor date
            # ------------------------------------------------------------------
            df_pay = discount_curve.df(payment_dt)

            if start_accrued_dt <= value_dt:
                # ----------------------------------------------------------------
                # Current (in-progress) period:
                #   - Equity price at period start is already fixed = current_price
                #     (or strike if no current price supplied and value_dt ==
                #     effective_dt).
                #   - We only need to forward from value_dt to end_accrued_dt.
                #   - Period-start notional is known: S_reset * quantity.
                # ----------------------------------------------------------------
                period_start_notional = self.current_price * self.quantity

                df_idx_end = index_curve.df(end_accrued_dt)
                df_div_end = dividend_curve.df(end_accrued_dt)
                df_idx_now = index_curve.df(value_dt)
                df_div_now = dividend_curve.df(value_dt)

                # Forward equity price from value_dt to end of period
                # F = S_now * (df_idx_now / df_idx_end) * (df_div_now / df_div_end)
                eq_growth = (df_idx_now / df_idx_end) * (df_div_now / df_div_end)
                fwd_end_notional = self.current_price * self.quantity * eq_growth

                payment_amount = fwd_end_notional - period_start_notional

                # Convenience rates for reporting (from value_dt to period end)
                index_alpha = self.year_fracs[i_pmnt]  # approximation; period may be partial
                fwd_rate    = (df_idx_now / df_idx_end - 1.0) / index_alpha if index_alpha > 0 else 0.0
                div_fwd_rate = (df_div_now / df_div_end - 1.0) / index_alpha if index_alpha > 0 else 0.0
                eq_fwd_rate  = (eq_growth - 1.0) / index_alpha if index_alpha > 0 else 0.0

            else:
                # ----------------------------------------------------------------
                # Future period: both t_{i-1} and t_i are beyond value_dt.
                #   Each such period is an independent forward contract.
                #
                #   PV_i = Q * S_0
                #          * [df_idx(0,t_{i-1}) * df_div(0,t_{i-1})]   <- growth to reset
                #          * df_disc(0,t_i)                              <- discount payment
                #          * [1/(df_idx(0,t_i)*df_div(0,t_i)) - 1]      <- net equity return
                #
                #   which equals:
                #     Q * S_0 * (df_idx_start * df_div_start / (df_idx_end * df_div_end) - 1)
                #             * df_disc(0,t_i)                           [after cancellation]
                # ----------------------------------------------------------------
                df_idx_start = index_curve.df(start_accrued_dt)
                df_idx_end   = index_curve.df(end_accrued_dt)
                df_div_start = dividend_curve.df(start_accrued_dt)
                df_div_end   = dividend_curve.df(end_accrued_dt)

                # Forward equity growth over the reset period
                eq_growth = (df_idx_start / df_idx_end) * (df_div_start / df_div_end)

                # Notional at period start (in expectation, funded from today)
                period_start_notional = self.current_price * self.quantity * df_idx_start * df_div_start

                # Expected notional at period end
                fwd_end_notional = period_start_notional * eq_growth  # = S_0 * Q * df_idx_start^2...
                # Rewrite cleanly:
                #   payment_amount (undiscounted, in forward measure) =
                #     S_0 * Q * (eq_growth_cumulative_to_end - eq_growth_cumulative_to_start)
                # Simplification: the forward payoff from the receiver's perspective is:
                #   S_{t_i} - S_{t_{i-1}}  (both unknown today)
                # Its PV = S_0 * Q * (df_idx_start*df_div_start/(df_idx_end*df_div_end) - 1) * df_pay/df_value
                #        = (fwd_price_end - fwd_price_start) * df_pay / df_value
                # where fwd_price_k = S_0 * df_idx(0,k) * df_div(0,k) / df_value [already in spot measure]

                fwd_price_start = self.current_price * (df_idx_start / df_value) * (df_div_start / df_value)
                fwd_price_end   = self.current_price * (df_idx_end   / df_value) * (df_div_end   / df_value)

                # *** Wait — cleaner canonical form (see note below) ***
                # PV_i = Q * [F(0,t_{i-1}) * df(0,t_i) - F(0,t_i) * df(0,t_i)]  ... doesn't simplify
                # Correct PV (standard result, no approximation):
                #   PV_i = Q * S_0 * df_div_start * df_idx_start   (long equity fwd to t_{i-1}, settled at t_i)
                #        - Q * S_0 * df_div_end   * df_idx_end      (short funded position at t_i)
                # Both already relative to curve anchor = value_dt base, so divide by df_value^2 ...
                # Actually if curves are anchored at value_dt, df(value_dt)=1 and this is exact:
                period_start_notional = self.current_price * self.quantity  # S_0 * Q, known today
                payment_amount = period_start_notional * (
                    (df_idx_start * df_div_start) / (df_idx_end * df_div_end) - 1.0
                )

                index_alpha  = self.year_fracs[i_pmnt]

                if index_alpha > 0:
                    fwd_rate = (df_idx_start / df_idx_end - 1.0) / index_alpha
                else:
                    fwd_rate = 0.0

                if index_alpha > 0:
                    div_fwd_rate = (df_div_start / df_div_end - 1.0) / index_alpha
                else:
                    div_fwd_rate = 0.0

                if index_alpha > 0:
                    eq_fwd_rate = ((df_idx_start * df_div_start) / (df_idx_end * df_div_end) - 1.0) / index_alpha
                else:
                    eq_fwd_rate = 0.0

            # Discount the payment — curves anchored at value_dt so df_value cancels
            df_payment = df_pay / df_value
            payment_pv = payment_amount * df_payment

            leg_pv += payment_pv

            self.fwd_rates.append(fwd_rate)
            self.div_fwd_rates.append(div_fwd_rate)
            self.eq_fwd_rates.append(eq_fwd_rate)
            self.last_notionals.append(period_start_notional)
            self.payment_amounts.append(payment_amount)
            self.payment_dfs.append(df_payment)
            self.payment_pvs.append(payment_pv)
            self.cumulative_pvs.append(leg_pv)

        if self.leg_type == SwapTypes.PAY:
            leg_pv *= -1.0

        return leg_pv

    ####################################################################################
    # THIS IS FROM A PULL REQUEST - I BELIEVE IT TO BE INCORRECT - IT ALREADY HAD A BUG
    ####################################################################################

    def value_bug(
        self,
        value_dt: Date,
        discount_curve: DiscountCurve,
        index_curve: DiscountCurve,
        dividend_curve: DiscountCurve = None,
        current_price: float = None
    ):
        """Value the equity leg with payments from an equity price, quantity,
        an index curve and an [optional] dividend curve. Discounting is based
        on a supplied discount curve as of the valuation date supplied.
        """

        if index_curve is None:
            index_curve = discount_curve

        if discount_curve is None:
            raise FinError("Discount curve not provided!")

        # Assume a naive dividend curve if nothing provided
        if dividend_curve is None:
            dividend_curve = FlatDiscountCurve(value_dt, 0)

        if discount_curve.value_dt != value_dt:
            raise FinError(
                "Discount Curve valuation date not same as value date"
            )

        # Current price can't be different from strike at effective date
        if current_price is not None:
            self.current_price = current_price
        else:
            self.current_price = self.strike

        self.fwd_rates = []
        self.div_fwd_rates = []
        self.eq_fwd_rates = []
        self.last_notionals = []
        self.payment_amounts = []
        self.payment_dfs = []
        self.payment_pvs = []
        self.cumulative_pvs = []

        df_value = discount_curve.df(value_dt)
        leg_pv = 0.0
        eq_term_rate = 0.0
        last_notional = self.notional # self.current_price * self.quantity
        next_notional = last_notional
        num_payments = len(self.payment_dts)

        index_day_counter = DayCount(self.accrual_dc_type)

        for i_pmnt in range(0, num_payments):

            payment_dt = self.payment_dts[i_pmnt]

            if payment_dt > value_dt:

                start_accrued_dt = self.start_accd_dts[i_pmnt]
                end_accrued_dt = self.end_accd_dts[i_pmnt]

                # CHANGED !!!! I BELIEVE I HAVE FIXED A BUG HERE FROM PR
                # Do not ask curves for dates before their valuation date.
                curve_start_dt = value_dt
                if start_accrued_dt > value_dt:
                    curve_start_dt = start_accrued_dt

                index_alpha = index_day_counter.year_frac(
                    curve_start_dt, end_accrued_dt
                )[0]

                df_start = index_curve.df(curve_start_dt)
                df_end = index_curve.df(end_accrued_dt)
                fwd_rate = (df_start / df_end - 1.0) / index_alpha

                div_start = dividend_curve.df(curve_start_dt)
                div_end = dividend_curve.df(end_accrued_dt)
                div_fwd_rate = (div_start / div_end - 1.0) / index_alpha

                # Equity discount derived from index and div curves
                eq_fwd_rate = (
                    (df_start / df_end) * (div_start / div_end) - 1
                ) / index_alpha

                self.fwd_rates.append(fwd_rate)
                self.div_fwd_rates.append(div_fwd_rate)
                self.eq_fwd_rates.append(eq_fwd_rate)

                # Iterative update of the term rate
                eq_term_rate = (1 + eq_fwd_rate * self.year_fracs[i_pmnt]) * (
                    1 + eq_term_rate
                ) - 1

                next_price = self.current_price * (1 + eq_term_rate)
                next_notional = next_price * self.quantity
                payment_amount = next_notional - last_notional

                df_payment = discount_curve.df(payment_dt) / df_value
                payment_pv = payment_amount * df_payment
                leg_pv += payment_pv

                self.last_notionals.append(last_notional)
                self.payment_amounts.append(payment_amount)
                self.payment_dfs.append(df_payment)
                self.payment_pvs.append(payment_pv)
                self.cumulative_pvs.append(leg_pv)

            else:

                self.fwd_rates.append(0.0)
                self.div_fwd_rates.append(0.0)
                self.eq_fwd_rates.append(0.0)
                self.last_notionals.append(self.notional)
                self.payment_amounts.append(0.0)
                self.payment_dfs.append(0.0)
                self.payment_pvs.append(0.0)
                self.cumulative_pvs.append(leg_pv)

            last_notional = next_notional

        if self.leg_type == SwapTypes.PAY:
            leg_pv = leg_pv * (-1.0)

        return leg_pv

    ##########################################################################

    def print_payment_amounts(self):
        """Prints the payment dates, accrual factors, discount factors,
        cash amounts, their present value and their cumulative PV using the
        last valuation performed."""

        print("START DATE:", self.effective_dt)
        print("MATURITY DATE:", self.maturity_dt)
        print("FREQUENCY:", str(self.freq_type))
        print("DAY COUNT:", str(self.accrual_dc_type))

        if len(self.payment_dts) == 0:
            print("Payments Dates not calculated.")
            return

        header = [
            "PAY_NUM",
            "PAY_dt",
            "ACCR_START",
            "ACCR_END",
            "DAYS",
            "YEARFRAC",
        ]

        rows = []
        num_flows = len(self.payment_dts)
        for i_flow in range(0, num_flows):
            rows.append(
                [
                    i_flow + 1,
                    self.payment_dts[i_flow],
                    self.start_accd_dts[i_flow],
                    self.end_accd_dts[i_flow],
                    self.accrued_days[i_flow],
                    round(self.year_fracs[i_flow], 4),
                ]
            )

        table = format_table(header, rows)
        print("\nPAYMENTS SCHEDULE:")
        print(table)

    ###########################################################################

    def print_valuation(self):
        """Prints the valuation dates, accrual factors, discount factors,
        cash amounts, their present value and their cumulative PV using the
        last valuation performed."""

        self.print_payment_amounts()

        if len(self.payment_amounts) == 0:
            print("Payments not calculated.")
            return

        header = [
            "PAY_NUM",
            "PAY_dt",
            "NOTIONAL",
            "FWD_RATE",
            "PMNT",
            "DF",
            "PV",
            "CUM_PV",
        ]

        rows = []
        num_flows = len(self.payment_dts)
        for i_flow in range(0, num_flows):
            rows.append(
                [
                    i_flow + 1,
                    self.payment_dts[i_flow],
                    round(self.last_notionals[i_flow], 0),
                    round(self.eq_fwd_rates[i_flow] * 100.0, 4),
                    round(self.payment_amounts[i_flow], 2),
                    round(self.payment_dfs[i_flow], 4),
                    round(self.payment_pvs[i_flow], 2),
                    round(self.cumulative_pvs[i_flow], 2),
                ]
            )

        table = format_table(header, rows)
        print("\nPAYMENTS VALUATION:")
        print(table)

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("EFFECTIVE DATE", self.effective_dt)
        s += label_to_string("MATURITY_DATE", self.maturity_dt)
        s += label_to_string("NOTIONAL", self.strike * self.quantity)
        s += label_to_string("SWAP TYPE", self.leg_type)
        s += label_to_string("RETURN TYPE", self.return_type)
        s += label_to_string("FREQUENCY", self.freq_type)
        s += label_to_string("DC_TYPE", self.accrual_dc_type)
        s += label_to_string("CALENDAR", self.cal_type)
        s += label_to_string("BUS DAY ADJUST", self.bd_type)
        s += label_to_string("DATE GEN TYPE", self.dg_type)
        return s


########################################################################################
