##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.math import ONE_MILLION
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
                 effective_date: Date,  # Contract starts or last equity reset
                 term_date_or_tenor: (Date, str),  # Date contract ends
                 leg_type: SwapTypes,
                 freq_type: FrequencyTypes,
                 day_count_type: DayCountTypes,
                 strike: float,  # Price at effective date
                 quantity: float = 1.0,  # Quantity at effective date
                 payment_lag: int = 0,
                 return_type: ReturnTypes = ReturnTypes.TOTAL_RETURN,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
                 end_of_month: bool = False):
        """ Create the equity leg of a swap contract giving the contract start
        date, its maturity, underlying strike price and quantity, payment
        frequency, day count convention, return type, and other details """

        check_argument_types(self.__init__, locals())

        if type(term_date_or_tenor) == Date:
            termination_date = term_date_or_tenor
        else:
            termination_date = effective_date.add_tenor(term_date_or_tenor)

        calendar = Calendar(calendar_type)
        self._maturity_date = calendar.adjust(termination_date,
                                              bus_day_adjust_type)

        if effective_date > self._maturity_date:
            raise FinError("Effective date after maturity date")

        if quantity < 0:
            # Long/Short is defined by leg_type
            raise FinError("Quantity must be non-negative")

        if return_type != ReturnTypes.TOTAL_RETURN:
            raise NotImplementedError("Return Type still not implemented")

        # To generate ISDA pmnt schedule properly we can't use these types
        if freq_type in (FrequencyTypes.CONTINUOUS, FrequencyTypes, FrequencyTypes.SIMPLE):
            raise FinError("Cannot generate payment schedule for this frequency!")

        self._effective_date = effective_date
        self._termination_date = termination_date
        self._leg_type = leg_type
        self._freq_type = freq_type
        self._payment_lag = payment_lag
        self._strike = strike
        self._quantity = quantity
        self._notional = strike * quantity
        self._return_type = return_type

        self._day_count_type = day_count_type
        self._calendar_type = calendar_type
        self._bus_day_adjust_type = bus_day_adjust_type
        self._date_gen_rule_type = date_gen_rule_type
        self._end_of_month = end_of_month

        self._startAccruedDates = []
        self._endAccruedDates = []
        self._payment_dates = []
        self._payments = []
        self._year_fracs = []
        self._accrued_days = []
        self._rates = []

        self.generate_payment_dates()

###############################################################################

    def generate_payment_dates(self):
        """ Generate the Equity leg payment dates and accrual factors. Similar
        to swap float leg, payment values can't be generated, as we do not have
        index curve, dividend curve and equity price. """

        schedule = Schedule(self._effective_date,
                            self._maturity_date,
                            self._freq_type,
                            self._calendar_type,
                            self._bus_day_adjust_type,
                            self._date_gen_rule_type,
                            end_of_month=self._end_of_month)

        scheduleDates = schedule._adjusted_dates

        if len(scheduleDates) < 2:
            raise FinError("Schedule has none or only one date")

        self._startAccruedDates = []
        self._endAccruedDates = []
        self._payment_dates = []
        self._year_fracs = []
        self._accrued_days = []

        prev_dt = scheduleDates[0]

        day_counter = DayCount(self._day_count_type)
        calendar = Calendar(self._calendar_type)

        # All of the lists end up with the same length
        for next_dt in scheduleDates[1:]:

            self._startAccruedDates.append(prev_dt)
            self._endAccruedDates.append(next_dt)

            if self._payment_lag == 0:
                payment_date = next_dt
            else:
                payment_date = calendar.add_business_days(next_dt,
                                                          self._payment_lag)

            self._payment_dates.append(payment_date)

            (year_frac, num, _) = day_counter.year_frac(prev_dt,
                                                        next_dt)

            self._year_fracs.append(year_frac)
            self._accrued_days.append(num)

            prev_dt = next_dt

###############################################################################

    def value(self,
              valuation_date: Date,
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

        if discount_curve._valuation_date != valuation_date:
            raise FinError(
                "Discount Curve valuation date not same as option value date")

        if index_curve is None:
            index_curve = discount_curve

        # Assume a naive dividend curve if nothing provided
        if dividend_curve is None:
            dividend_curve = DiscountCurveFlat(valuation_date, 0)

        # Current price can't be different than strike at effective date
        if current_price is not None:
            self._current_price = current_price
        else:
            self._current_price = self._strike

        self._fwd_rates = []
        self._div_fwd_rates = []
        self._eq_fwd_rates = []
        self._last_notionals = []
        self._payments = []
        self._paymentDfs = []
        self._paymentPVs = []
        self._cumulativePVs = []

        dfValue = discount_curve.df(valuation_date)
        legPV, eq_term_rate = 0.0, 0.0
        lastNotional = self._notional
        numPayments = len(self._payment_dates)

        index_basis = index_curve._day_count_type
        index_day_counter = DayCount(index_basis)

        for iPmnt in range(0, numPayments):

            pmntDate = self._payment_dates[iPmnt]

            if pmntDate > valuation_date:

                startAccruedDt = self._startAccruedDates[iPmnt]
                endAccruedDt = self._endAccruedDates[iPmnt]
                index_alpha = index_day_counter.year_frac(startAccruedDt,
                                                          endAccruedDt)[0]

                dfStart = index_curve.df(startAccruedDt)
                dfEnd = index_curve.df(endAccruedDt)
                fwd_rate = (dfStart / dfEnd - 1.0) / index_alpha

                divStart = dividend_curve.df(startAccruedDt)
                divEnd = dividend_curve.df(endAccruedDt)
                div_fwd_rate = (divStart / divEnd - 1.0) / index_alpha

                # Equity discount derived from index and div curves
                eq_fwd_rate = ((dfStart / dfEnd) * (divStart / divEnd) - 1.0) / index_alpha

                self._fwd_rates.append(fwd_rate)
                self._div_fwd_rates.append(div_fwd_rate)
                self._eq_fwd_rates.append(eq_fwd_rate)

                # Iterative update of the term rate
                eq_term_rate = (1 + eq_fwd_rate * self._year_fracs[iPmnt]) * (1 + eq_term_rate)  - 1

                nextPrice = self._current_price * (1 + eq_term_rate)
                nextNotional = nextPrice * self._quantity
                pmntAmount = nextNotional - lastNotional

                dfPmnt = discount_curve.df(pmntDate) / dfValue
                pmntPV = pmntAmount * dfPmnt
                legPV += pmntPV

                self._last_notionals.append(lastNotional)
                self._payments.append(pmntAmount)
                self._paymentDfs.append(dfPmnt)
                self._paymentPVs.append(pmntPV)
                self._cumulativePVs.append(legPV)

            else:

                self._fwd_rates.append(0.0)
                self._div_fwd_rates.append(0.0)
                self._eq_fwd_rates.append(0.0)
                self._last_notionals.append(self._notional)
                self._payments.append(0.0)
                self._paymentDfs.append(0.0)
                self._paymentPVs.append(0.0)
                self._cumulativePVs.append(legPV)

            lastNotional = nextNotional

        if self._leg_type == SwapTypes.PAY:
            legPV = legPV * (-1.0)

        return legPV

##########################################################################

    def print_payments(self):
        """ Prints the payment dates, accrual factors, discount factors,
        cash amounts, their present value and their cumulative PV using the
        last valuation performed. """

        print("START DATE:", self._effective_date)
        print("MATURITY DATE:", self._maturity_date)
        print("FREQUENCY:", str(self._freq_type))
        print("DAY COUNT:", str(self._day_count_type))

        if len(self._payment_dates) == 0:
            print("Payments Dates not calculated.")
            return

        header = ["PAY_NUM", "PAY_DATE", "ACCR_START", "ACCR_END", "DAYS",
                  "YEARFRAC"]

        rows = []
        num_flows = len(self._payment_dates)
        for iFlow in range(0, num_flows):
            rows.append([
                iFlow + 1,
                self._payment_dates[iFlow],
                self._startAccruedDates[iFlow],
                self._endAccruedDates[iFlow],
                self._accrued_days[iFlow],
                round(self._year_fracs[iFlow], 4),
            ])

        table = format_table(header, rows)
        print("\nPAYMENTS SCHEDULE:")
        print(table)

###############################################################################

    def print_valuation(self):
        """ Prints the valuation dates, accrual factors, discount factors,
        cash amounts, their present value and their cumulative PV using the
        last valuation performed. """

        self.print_payments()

        if len(self._payments) == 0:
            print("Payments not calculated.")
            return

        header = ["PAY_NUM", "PAY_DATE",  "NOTIONAL",
                  "FWD_RATE", "PMNT", "DF",
                  "PV", "CUM_PV"]

        rows = []
        num_flows = len(self._payment_dates)
        for iFlow in range(0, num_flows):
            rows.append([
                iFlow + 1,
                self._payment_dates[iFlow],
                round(self._last_notionals[iFlow], 0),
                round(self._eq_fwd_rates[iFlow] * 100.0, 4),
                round(self._payments[iFlow], 2),
                round(self._paymentDfs[iFlow], 4),
                round(self._paymentPVs[iFlow], 2),
                round(self._cumulativePVs[iFlow], 2),
            ])

        table = format_table(header, rows)
        print("\nPAYMENTS VALUATION:")
        print(table)

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EFFECTIVE DATE", self._effective_date)
        s += label_to_string("MATURITY DATE", self._maturity_date)
        s += label_to_string("NOTIONAL", self._strike * self._quantity)
        s += label_to_string("SWAP TYPE", self._leg_type)
        s += label_to_string("RETURN TYPE", self._return_type)
        s += label_to_string("FREQUENCY", self._freq_type)
        s += label_to_string("DAY COUNT", self._day_count_type)
        s += label_to_string("CALENDAR", self._calendar_type)
        s += label_to_string("BUS DAY ADJUST", self._bus_day_adjust_type)
        s += label_to_string("DATE GEN TYPE", self._date_gen_rule_type)
        return s


###############################################################################
