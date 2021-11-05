##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.global_vars import gSmall
from ...utils.day_count import DayCountTypes
from ...utils.frequency import FrequencyTypes, annual_frequency
from ...utils.calendar import CalendarTypes, DateGenRuleTypes
from ...utils.calendar import Calendar, BusDayAdjustTypes
from ...utils.helpers import check_argument_types, label_to_string
from ...utils.math import ONE_MILLION
from ...utils.global_types import SwapTypes
from ...market.curves.discount_curve import DiscountCurve

from .swap_fixed_leg import SwapFixedLeg
from .swap_float_leg import SwapFloatLeg


##########################################################################


class IborSwap:
    """ Class for managing a standard Fixed vs IBOR swap. This is a contract
    in which a fixed payment leg is exchanged for a series of floating rates
    payments linked to some IBOR index rate. There is no exchange of principal.
    The contract is entered into at zero initial cost. The contract lasts from
    a start date to a specified maturity date.

    The floating rate is not known fully until the end of the preceding payment
    period. It is set in advance and paid in arrears.

    The value of the contract is the NPV of the two coupon streams. Discounting
    is done on a supplied discount curve which is separate from the curve from
    which the implied index rates are extracted. """

    def __init__(self,
                 effective_date: Date,  # Date interest starts to accrue
                 termination_date_or_tenor: (Date, str),  # Date contract ends
                 fixed_leg_type: SwapTypes,
                 fixed_coupon: float,  # Fixed coupon (annualised)
                 fixed_freq_type: FrequencyTypes,
                 fixed_day_count_type: DayCountTypes,
                 notional: float = ONE_MILLION,
                 float_spread: float = 0.0,
                 float_freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 float_day_count_type: DayCountTypes = DayCountTypes.THIRTY_E_360,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):
        """ Create an interest rate swap contract giving the contract start
        date, its maturity, fixed coupon, fixed leg frequency, fixed leg day
        count convention and notional. The floating leg parameters have default
        values that can be overwritten if needed. The start date is contractual
        and is the same as the settlement date for a new swap. It is the date
        on which interest starts to accrue. The end of the contract is the
        termination date. This is not adjusted for business days. The adjusted
        termination date is called the maturity date. This is calculated. """

        check_argument_types(self.__init__, locals())

        if type(termination_date_or_tenor) == Date:
            self._termination_date = termination_date_or_tenor
        else:
            self._termination_date = effective_date.add_tenor(
                termination_date_or_tenor)

        calendar = Calendar(calendar_type)
        self._maturity_date = calendar.adjust(self._termination_date,
                                              bus_day_adjust_type)

        if effective_date > self._maturity_date:
            raise FinError("Start date after maturity date")

        self._effective_date = effective_date

        float_leg_type = SwapTypes.PAY
        if fixed_leg_type == SwapTypes.PAY:
            float_leg_type = SwapTypes.RECEIVE

        payment_lag = 0
        principal = 0.0

        self._fixed_leg = SwapFixedLeg(effective_date,
                                       self._termination_date,
                                       fixed_leg_type,
                                       fixed_coupon,
                                       fixed_freq_type,
                                       fixed_day_count_type,
                                       notional,
                                       principal,
                                       payment_lag,
                                       calendar_type,
                                       bus_day_adjust_type,
                                       date_gen_rule_type)

        self._float_leg = SwapFloatLeg(effective_date,
                                       self._termination_date,
                                       float_leg_type,
                                       float_spread,
                                       float_freq_type,
                                       float_day_count_type,
                                       notional,
                                       principal,
                                       payment_lag,
                                       calendar_type,
                                       bus_day_adjust_type,
                                       date_gen_rule_type)

    ###########################################################################

    def value(self,
              valuation_date: Date,
              discount_curve: DiscountCurve,
              index_curve: DiscountCurve = None,
              firstFixingRate=None):
        """ Value the interest rate swap on a value date given a single Ibor
        discount curve. """

        if index_curve is None:
            index_curve = discount_curve

        fixed_leg_value = self._fixed_leg.value(valuation_date,
                                                discount_curve)

        float_leg_value = self._float_leg.value(valuation_date,
                                                discount_curve,
                                                index_curve,
                                                firstFixingRate)

        value = fixed_leg_value + float_leg_value
        return value

    ###########################################################################

    def pv01(self, valuation_date, discount_curve):
        """ Calculate the value of 1 basis point coupon on the fixed leg. """

        pv = self._fixed_leg.value(valuation_date, discount_curve)
        pv01 = pv / self._fixed_leg._coupon / self._fixed_leg._notional
        # Needs to be positive even if it is a payer leg
        pv01 = np.abs(pv01)
        return pv01

    ###########################################################################

    def swap_rate(self,
                  valuation_date: Date,
                  discount_curve: DiscountCurve,
                  index_curve: DiscountCurve = None,
                  first_fixing: float = None):
        """ Calculate the fixed leg coupon that makes the swap worth zero.
        If the valuation date is before the swap payments start then this
        is the forward swap rate as it starts in the future. The swap rate
        is then a forward swap rate and so we use a forward discount
        factor. If the swap fixed leg has begun then we have a spot
        starting swap. The swap rate can also be calculated in a dual curve
        approach but in this case the first fixing on the floating leg is
        needed. """

        pv01 = self.pv01(valuation_date, discount_curve)

        if abs(pv01) < gSmall:
            raise FinError("PV01 is zero. Cannot compute swap rate.")

        if valuation_date < self._effective_date:
            df0 = discount_curve.df(self._effective_date)
        else:
            df0 = discount_curve.df(valuation_date)

        float_leg_pv = 0.0

        if index_curve is None:
            df_t = discount_curve.df(self._maturity_date)
            float_leg_pv = (df0 - df_t)
        else:
            float_leg_pv = self._float_leg.value(valuation_date,
                                                 discount_curve,
                                                 index_curve,
                                                 first_fixing)

            float_leg_pv /= self._fixed_leg._notional

        cpn = float_leg_pv / pv01
        return cpn

    ###########################################################################

    def cash_settled_pv01(self,
                          valuation_date,
                          flat_swap_rate,
                          frequency_type):
        """ Calculate the forward value of an annuity of a forward starting
        swap using a single flat discount rate equal to the swap rate. This is
        used in the pricing of a cash-settled swaption in the IborSwaption
        class. This method does not affect the standard valuation methods."""

        m = annual_frequency(frequency_type)

        if m == 0:
            raise FinError("Frequency cannot be zero.")

        """ The swap may have started in the past but we can only value
        payments that have occurred after the valuation date. """
        start_index = 0
        while self._fixed_leg._payment_dates[start_index] < valuation_date:
            start_index += 1

        """ If the swap has yet to settle then we do not include the
        start date of the swap as a coupon payment date. """
        if valuation_date <= self._effective_date:
            start_index = 1

        """ Now PV fixed leg flows. """
        flatPV01 = 0.0
        df = 1.0
        alpha = 1.0 / m

        for _ in self._fixed_leg._payment_dates[start_index:]:
            df = df / (1.0 + alpha * flat_swap_rate)
            flatPV01 += df * alpha

        return flatPV01

    ###########################################################################

    def print_fixed_leg_pv(self):
        """ Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows. """

        self._fixed_leg.print_valuation()

    ###########################################################################

    def print_float_leg_pv(self):
        """ Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows. """

        self._float_leg.print_valuation()

    ###########################################################################

    def print_flows(self):
        """ Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows. """

        self._fixed_leg.print_payments()
        self._float_leg.print_payments()

    ###########################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += self._fixed_leg.__repr__()
        s += "\n"
        s += self._float_leg.__repr__()
        return s

    ###########################################################################

    def _print(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        print(self)

###############################################################################
