##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.calendar import CalendarTypes,  DateGenRuleTypes
from ...utils.calendar import Calendar, BusDayAdjustTypes
from ...utils.schedule import Schedule
from ...utils.helpers import format_table, label_to_string, check_argument_types
from ...utils.global_types import SwapTypes, ReturnTypes
from ...market.curves.discount_curve import DiscountCurve
from ...market.curves.discount_curve_flat import DiscountCurveFlat

##########################################################################


class EquitySwapLeg:
    """ Class for managing the equity leg of an equity swap. An equity leg is
    a leg with a sequence of flows calculated according to an ISDA schedule
    and follows the economics of a collection of equity forward contracts.
    """

    def __init__(self,
                 effective_dt: Date,  # Contract starts or last equity reset
                 term_dt_or_tenor: (Date, str),  # Date contract ends
                 leg_type: SwapTypes,
                 freq_type: FrequencyTypes,
                 dc_type: DayCountTypes,
                 strike: float,  # Price at effective date
                 quantity: float = 1.0,  # Quantity at effective date
                 payment_lag: int = 0,
                 return_type: ReturnTypes = ReturnTypes.TOTAL_RETURN,
                 cal_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
                 end_of_month: bool = False):
        """ Create the equity leg of a swap contract giving the contract start
        date, its maturity, underlying strike price and quantity, payment
        frequency, day count convention, return type, and other details """

        check_argument_types(self.__init__, locals())

        if isinstance(term_dt_or_tenor, Date):
            termination_dt = term_dt_or_tenor
        else:
            termination_dt = effective_dt.add_tenor(term_dt_or_tenor)

        calendar = Calendar(cal_type)
        self._maturity_dt = calendar.adjust(termination_dt,
                                            bd_type)

        if effective_dt > self._maturity_dt:
            raise FinError("Effective date after maturity date")

        if quantity < 0:
            # Long/Short is defined by leg_type
            raise FinError("Quantity must be non-negative")

        if return_type != ReturnTypes.TOTAL_RETURN:
            raise NotImplementedError("Return Type still not implemented")

        # To generate ISDA pmnt schedule properly we can't use these types
        if freq_type in (FrequencyTypes.CONTINUOUS, FrequencyTypes,
                         FrequencyTypes.SIMPLE):
            raise FinError("Cannot generate payment schedule for this frequency!")

        self._effective_dt = effective_dt
        self._termination_dt = termination_dt
        self._leg_type = leg_type
        self._freq_type = freq_type
        self._payment_lag = payment_lag
        self._strike = strike
        self._quantity = quantity
        self._notional = strike * quantity
        self._return_type = return_type

        self._dc_type = dc_type
        self._cal_type = cal_type
        self._bd_type = bd_type
        self._dg_type = dg_type
        self._end_of_month = end_of_month

        self._start_accd_dts = []
        self._end_accd_dts = []
        self._payment_dts = []
        self._pmnts = []
        self._year_fracs = []
        self._accrued_days = []
        self._rates = []

        self.generate_payment_dts()

###############################################################################

    def generate_payment_dts(self):
        """ Generate the Equity leg payment dates and accrual factors. Similar
        to swap float leg, payment values can't be generated, as we do not have
        index curve, dividend curve and equity price. """

        schedule = Schedule(self._effective_dt,
                            self._maturity_dt,
                            self._freq_type,
                            self._cal_type,
                            self._bd_type,
                            self._dg_type,
                            end_of_month=self._end_of_month)

        schedule_dts = schedule._adjusted_dts

        if len(schedule_dts) < 2:
            raise FinError("Schedule has none or only one date")

        self._start_accd_dts = []
        self._end_accd_dts = []
        self._payment_dts = []
        self._year_fracs = []
        self._accrued_days = []

        prev_dt = schedule_dts[0]

        day_counter = DayCount(self._dc_type)
        calendar = Calendar(self._cal_type)

        # All of the lists end up with the same length
        for next_dt in schedule_dts[1:]:

            self._start_accd_dts.append(prev_dt)
            self._end_accd_dts.append(next_dt)

            if self._payment_lag == 0:
                payment_dt = next_dt
            else:
                payment_dt = calendar.add_business_days(next_dt,
                                                        self._payment_lag)

            self._payment_dts.append(payment_dt)

            (year_frac, num, _) = day_counter.year_frac(prev_dt,
                                                        next_dt)

            self._year_fracs.append(year_frac)
            self._accrued_days.append(num)

            prev_dt = next_dt

###############################################################################

    def value(self,
              value_dt: Date,
              discount_curve: DiscountCurve,
              index_curve: DiscountCurve,
              dividend_curve: DiscountCurve = None,
              current_price: float = None):
        """ Value the equity leg with payments from an equity price, quantity,
        an index curve and an [optional] dividend curve. Discounting is based
        on a supplied discount curve as of the valuation date supplied.
        """

        if discount_curve is None:
            raise FinError("Discount curve not provided!")

        if discount_curve.value_dt != value_dt:
            raise FinError(
                "Discount Curve valuation date not same as option value date")

        if index_curve is None:
            index_curve = discount_curve

        # Assume a naive dividend curve if nothing provided
        if dividend_curve is None:
            dividend_curve = DiscountCurveFlat(value_dt, 0)

        # Current price can't be different than strike at effective date
        if current_price is not None:
            self._current_price = current_price
        else:
            self._current_price = self._strike

        self._fwd_rates = []
        self._div_fwd_rates = []
        self._eq_fwd_rates = []
        self._last_notionals = []
        self._pmnts = []
        self._pmnt_dfs = []
        self._pmnt_pvs = []
        self._cumulative_pvs = []

        df_value = discount_curve.df(value_dt)
        leg_pv, eq_term_rate = 0.0, 0.0
        last_notional = self._notional
        num_payments = len(self._payment_dts)

        index_basis = index_curve._dc_type
        index_day_counter = DayCount(index_basis)

        for i_pmnt in range(0, num_payments):

            pmnt_dt = self._payment_dts[i_pmnt]

            if pmnt_dt > value_dt:

                start_accrued_dt = self._start_accd_dts[i_pmnt]
                end_accrued_dt = self._end_accd_dts[i_pmnt]
                index_alpha = index_day_counter.year_frac(start_accrued_dt,
                                                          end_accrued_dt)[0]

                df_start = index_curve.df(start_accrued_dt)
                df_end = index_curve.df(end_accrued_dt)
                fwd_rate = (df_start / df_end - 1.0) / index_alpha

                div_start = dividend_curve.df(start_accrued_dt)
                div_end = dividend_curve.df(end_accrued_dt)
                div_fwd_rate = (div_start / div_end - 1.0) / index_alpha

                # Equity discount derived from index and div curves
                eq_fwd_rate = ((df_start / df_end) * (div_start / div_end) - 1) / index_alpha

                self._fwd_rates.append(fwd_rate)
                self._div_fwd_rates.append(div_fwd_rate)
                self._eq_fwd_rates.append(eq_fwd_rate)

                # Iterative update of the term rate
                eq_term_rate = (1 + eq_fwd_rate * self._year_fracs[i_pmnt]) * (1 + eq_term_rate)  - 1

                next_price = self._current_price * (1 + eq_term_rate)
                next_notional = next_price * self._quantity
                pmntAmount = next_notional - last_notional

                df_pmnt = discount_curve.df(pmnt_dt) / df_value
                pmnt_pv = pmntAmount * df_pmnt
                leg_pv += pmnt_pv

                self._last_notionals.append(last_notional)
                self._pmnts.append(pmntAmount)
                self._pmnt_dfs.append(df_pmnt)
                self._pmnt_pvs.append(pmnt_pv)
                self._cumulative_pvs.append(leg_pv)

            else:

                self._fwd_rates.append(0.0)
                self._div_fwd_rates.append(0.0)
                self._eq_fwd_rates.append(0.0)
                self._last_notionals.append(self._notional)
                self._pmnts.append(0.0)
                self._pmnt_dfs.append(0.0)
                self._pmnt_pvs.append(0.0)
                self._cumulative_pvs.append(leg_pv)

            last_notional = next_notional

        if self._leg_type == SwapTypes.PAY:
            leg_pv = leg_pv * (-1.0)

        return leg_pv

##########################################################################

    def print_pmnts(self):
        """ Prints the payment dates, accrual factors, discount factors,
        cash amounts, their present value and their cumulative PV using the
        last valuation performed. """

        print("START DATE:", self._effective_dt)
        print("MATURITY DATE:", self._maturity_dt)
        print("FREQUENCY:", str(self._freq_type))
        print("DAY COUNT:", str(self._dc_type))

        if len(self._payment_dts) == 0:
            print("Payments Dates not calculated.")
            return

        header = ["PAY_NUM", "PAY_dt", "ACCR_START", "ACCR_END", "DAYS",
                  "YEARFRAC"]

        rows = []
        num_flows = len(self._payment_dts)
        for i_flow in range(0, num_flows):
            rows.append([
                i_flow + 1,
                self._payment_dts[i_flow],
                self._start_accd_dts[i_flow],
                self._end_accd_dts[i_flow],
                self._accrued_days[i_flow],
                round(self._year_fracs[i_flow], 4),
            ])

        table = format_table(header, rows)
        print("\nPAYMENTS SCHEDULE:")
        print(table)

###############################################################################

    def print_valuation(self):
        """ Prints the valuation dates, accrual factors, discount factors,
        cash amounts, their present value and their cumulative PV using the
        last valuation performed. """

        self.print_pmnts()

        if len(self._pmnts) == 0:
            print("Payments not calculated.")
            return

        header = ["PAY_NUM", "PAY_dt",  "NOTIONAL",
                  "FWD_RATE", "PMNT", "DF",
                  "PV", "CUM_PV"]

        rows = []
        num_flows = len(self._payment_dts)
        for i_flow in range(0, num_flows):
            rows.append([
                i_flow + 1,
                self._payment_dts[i_flow],
                round(self._last_notionals[i_flow], 0),
                round(self._eq_fwd_rates[i_flow] * 100.0, 4),
                round(self._pmnts[i_flow], 2),
                round(self._pmnt_dfs[i_flow], 4),
                round(self._pmnt_pvs[i_flow], 2),
                round(self._cumulative_pvs[i_flow], 2),
            ])

        table = format_table(header, rows)
        print("\nPAYMENTS VALUATION:")
        print(table)

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EFFECTIVE DATE", self._effective_dt)
        s += label_to_string("MATURITY DATE", self._maturity_dt)
        s += label_to_string("NOTIONAL", self._strike * self._quantity)
        s += label_to_string("SWAP TYPE", self._leg_type)
        s += label_to_string("RETURN TYPE", self._return_type)
        s += label_to_string("FREQUENCY", self._freq_type)
        s += label_to_string("DAY COUNT", self._dc_type)
        s += label_to_string("CALENDAR", self._cal_type)
        s += label_to_string("BUS DAY ADJUST", self._bd_type)
        s += label_to_string("DATE GEN TYPE", self._dg_type)
        return s


###############################################################################
