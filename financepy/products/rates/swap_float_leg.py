###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################
import numpy as np
import pandas as pd

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.math import ONE_MILLION
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
from ...utils.global_types import SwapTypes
from ...market.curves.discount_curve import DiscountCurve

###############################################################################


class SwapFloatLeg:
    """Class for managing the floating leg of a swap. A float leg consists of
    a sequence of flows calculated according to an ISDA schedule and with a
    coupon determined by an index curve which changes over life of the swap."""

    def __init__(
        self,
        effective_dt: Date,  # Date interest starts to accrue
        end_dt: (Date, str),  # Date contract ends
        leg_type: SwapTypes,
        spread: float,
        freq_type: FrequencyTypes,
        dc_type: DayCountTypes,
        notional: float = ONE_MILLION,
        principal: float = 0.0,
        payment_lag: int = 0,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
        end_of_month: bool = False,
    ):
        """Create the fixed leg of a swap contract giving the contract start
        date, its maturity, fixed coupon, fixed leg frequency, fixed leg day
        count convention and notional."""

        check_argument_types(self.__init__, locals())

        if type(end_dt) is Date:
            self.termination_dt = end_dt
        else:
            self.termination_dt = effective_dt.add_tenor(end_dt)

        calendar = Calendar(cal_type)

        self.maturity_dt = calendar.adjust(self.termination_dt, bd_type)

        if effective_dt > self.maturity_dt:
            raise FinError("Start date after maturity date")

        self.effective_dt = effective_dt
        self.end_dt = end_dt
        self.leg_type = leg_type
        self.freq_type = freq_type
        self.payment_lag = payment_lag
        self.principal = 0.0
        self.notional = notional
        self.notional_array = []
        self.spread = spread

        self.dc_type = dc_type
        self.cal_type = cal_type
        self.bd_type = bd_type
        self.dg_type = dg_type
        self.end_of_month = end_of_month

        self.start_accrued_dts = []
        self.end_accrued_dts = []
        self.payment_dts = []
        self.payments = []
        self.year_fracs = []
        self.accrued_days = []

        self.generate_payment_dts()

    ###########################################################################

    def generate_payment_dts(self):
        """Generate the floating leg payment dates and accrual factors. The
        coupons cannot be generated yet as we do not have the index curve."""

        schedule = Schedule(
            self.effective_dt,
            self.termination_dt,
            self.freq_type,
            self.cal_type,
            self.bd_type,
            self.dg_type,
            end_of_month=self.end_of_month,
        )

        schedule_dts = schedule.adjusted_dts

        if len(schedule_dts) < 2:
            raise FinError("Schedule has none or only one date")

        self.start_accrued_dts = []
        self.end_accrued_dts = []
        self.payment_dts = []
        self.year_fracs = []
        self.accrued_days = []

        prev_dt = schedule_dts[0]

        day_counter = DayCount(self.dc_type)
        calendar = Calendar(self.cal_type)

        # All of the lists end up with the same length
        for next_dt in schedule_dts[1:]:

            self.start_accrued_dts.append(prev_dt)
            self.end_accrued_dts.append(next_dt)

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
        value_dt: Date,  # This should be the settlement date
        discount_curve: DiscountCurve,
        index_curve: DiscountCurve,
        first_fixing_rate: float = None,
        pv_only=True,
    ):
        """Value the floating leg with payments from an index curve and
        discounting based on a supplied discount curve as of the valuation date
        supplied. For an existing swap, the user must enter the next fixing
        coupon."""

        if discount_curve is None:
            raise FinError("Discount curve is None")

        if index_curve is None:
            index_curve = discount_curve

        self.rates = []
        self.payments = []
        self.payment_dfs = []
        self.payment_pvs = []
        self.cumulative_pvs = []

        df_value = discount_curve.df(value_dt)
        leg_pv = 0.0
        num_payments = len(self.payment_dts)
        first_payment = False

        if not len(self.notional_array):
            self.notional_array = [self.notional] * num_payments

        index_basis = index_curve.dc_type
        index_day_counter = DayCount(index_basis)

        for i_pmnt in range(0, num_payments):

            payment_dt = self.payment_dts[i_pmnt]

            if payment_dt > value_dt:

                start_accrued_dt = self.start_accrued_dts[i_pmnt]
                end_accrued_dt = self.end_accrued_dts[i_pmnt]
                pay_alpha = self.year_fracs[i_pmnt]

                (index_alpha, num, _) = index_day_counter.year_frac(
                    start_accrued_dt, end_accrued_dt
                )

                if first_payment is False and first_fixing_rate is not None:

                    fwd_rate = first_fixing_rate
                    first_payment = True

                else:

                    df_start = index_curve.df(start_accrued_dt)
                    df_end = index_curve.df(end_accrued_dt)
                    fwd_rate = (df_start / df_end - 1.0) / index_alpha

                payment_amount = (
                    (fwd_rate + self.spread)
                    * pay_alpha
                    * self.notional_array[i_pmnt]
                )

                df_payment = discount_curve.df(payment_dt) / df_value
                payment_pv = payment_amount * df_payment
                leg_pv += payment_pv

                self.rates.append(fwd_rate)
                self.payments.append(payment_amount)
                self.payment_dfs.append(df_payment)
                self.payment_pvs.append(payment_pv)
                self.cumulative_pvs.append(leg_pv)

            else:

                self.rates.append(0.0)
                self.payments.append(0.0)
                self.payment_dfs.append(0.0)
                self.payment_pvs.append(0.0)
                self.cumulative_pvs.append(leg_pv)

        if payment_dt > value_dt:
            payment_pv = self.principal * df_payment * self.notional_array[-1]
            self.payment_pvs[-1] += payment_pv
            leg_pv += payment_pv
            self.cumulative_pvs[-1] = leg_pv

        if self.leg_type == SwapTypes.PAY:
            leg_pv = leg_pv * (-1.0)

        if pv_only:
            return leg_pv
        else:
            return leg_pv, self._cashflow_report_from_cached_values()

    ###########################################################################

    def _cashflow_report_from_cached_values(self):
        """After calling value(...) function, internal members store
        cashflow-by-cashflow values
        Return them in a dataframe

        Returns:
            pd.DataFrame: cashflow values and related data
        """

        leg_type_sign = -1 if self.leg_type == SwapTypes.PAY else 1
        df = pd.DataFrame()
        df["payment_date"] = self.payment_dts
        df["start_accrual_date"] = self.start_accrued_dts
        df["end_accrual_date"] = self.end_accrued_dts
        df["year_frac"] = self.year_fracs
        df["rate"] = self.rates
        df["payment"] = np.array(self.payments) * leg_type_sign
        df["payment_df"] = self.payment_dfs
        df["payment_pv"] = np.array(self.payment_pvs) * leg_type_sign
        df["leg"] = "FLOAT"

        return df

    ###########################################################################

    def print_payments(self):
        """Prints the fixed leg dates, accrual factors, discount factors,
        cash amounts, their present value and their cumulative PV using the
        last valuation performed."""

        print("START DATE:", self.effective_dt)
        print("MATURITY DATE:", self.maturity_dt)
        print("SPREAD (bp):", self.spread * 10000)
        print("FREQUENCY:", str(self.freq_type))
        print("DAY COUNT:", str(self.dc_type))

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
                    self.start_accrued_dts[i_flow],
                    self.end_accrued_dts[i_flow],
                    self.accrued_days[i_flow],
                    round(self.year_fracs[i_flow], 4),
                ]
            )

        table = format_table(header, rows)
        print("\nPAYMENTS SCHEDULE:")
        print(table)

    ###########################################################################

    def print_valuation(self):
        """Prints the fixed leg dates, accrual factors, discount factors,
        cash amounts, their present value and their cumulative PV using the
        last valuation performed."""

        print("START DATE:", self.effective_dt)
        print("MATURITY DATE:", self.maturity_dt)
        print("SPREAD (BPS):", self.spread * 10000)
        print("FREQUENCY:", str(self.freq_type))
        print("DAY COUNT:", str(self.dc_type))

        if len(self.payments) == 0:
            print("Payments not calculated.")
            return

        header = [
            "PAY_NUM",
            "PAY_dt",
            "NOTIONAL",
            "IBOR",
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
                    round(self.notional_array[i_flow], 0),
                    round(self.rates[i_flow] * 100.0, 4),
                    round(self.payments[i_flow], 2),
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
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("START DATE", self.effective_dt)
        s += label_to_string("TERMINATION DATE", self.termination_dt)
        s += label_to_string("MATURITY DATE", self.maturity_dt)
        s += label_to_string("NOTIONAL", self.notional)
        s += label_to_string("SWAP TYPE", self.leg_type)
        s += label_to_string("SPREAD (BPS)", self.spread * 10000)
        s += label_to_string("FREQUENCY", self.freq_type)
        s += label_to_string("DAY COUNT", self.dc_type)
        s += label_to_string("CALENDAR", self.cal_type)
        s += label_to_string("BUS DAY ADJUST", self.bd_type)
        s += label_to_string("DATE GEN TYPE", self.dg_type)
        return s

    ###########################################################################

    def _print(self):
        """Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond."""
        print(self)


###############################################################################
