##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.day_count import DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.calendar import CalendarTypes,  DateGenRuleTypes
from ...utils.calendar import Calendar, BusDayAdjustTypes
from ...utils.helpers import check_argument_types, label_to_string
from ...utils.math import ONE_MILLION
from ...utils.global_types import SwapTypes
from ...market.curves.discount_curve import DiscountCurve

from .swap_fixed_leg import SwapFixedLeg
from .swap_float_leg import SwapFloatLeg

###############################################################################

from enum import Enum


class FinCompoundingTypes(Enum):
    COMPOUNDED = 1
    OVERNIGHT_COMPOUNDED_ANNUAL_RATE = 2
    AVERAGED = 3
    AVERAGED_DAILY = 4


###############################################################################

class OIS:
    """ Class for managing overnight index rate swaps (OIS) and Fed Fund swaps. 
    This is a contract in which a fixed payment leg is exchanged for a payment
    which pays the rolled-up overnight index rate (OIR). There is no exchange
    of par. The contract is entered into at zero initial cost.

    NOTE: This class is almost identical to IborSwap but will possibly
    deviate as distinctions between the two become clear to me. If not they 
    will be converged (or inherited) to avoid duplication.

    The contract lasts from a start date to a specified maturity date.
    The fixed coupon is the OIS fixed rate for the corresponding tenor which is
    set at contract initiation.

    The floating rate is not known fully until the end of each payment period.
    It's calculated at the contract maturity and is based on daily observations
    of the overnight index rate which are compounded according to a specific
    convention. Hence the OIS floating rate is determined by the history of the
    OIS rates.

    In its simplest form, there is just one fixed rate payment and one floating
    rate payment at contract maturity. However when the contract becomes longer
    than one year the floating and fixed payments become periodic, usually with
    annual exchanges of cash.

    The value of the contract is the NPV of the two coupon streams. Discounting
    is done on the OIS curve which is itself implied by the term structure of
    market OIS rates. """

    def __init__(self,
                 effective_date: Date,  # Date interest starts to accrue
                 termination_date_or_tenor: (Date, str),  # Date contract ends
                 fixed_leg_type: SwapTypes,
                 fixed_coupon: float,  # Fixed coupon (annualised)
                 fixed_freq_type: FrequencyTypes,
                 fixed_day_count_type: DayCountTypes,
                 notional: float = ONE_MILLION,
                 payment_lag: int = 0,  # Number of days after period payment occurs
                 float_spread: float = 0.0,
                 float_freq_type: FrequencyTypes = FrequencyTypes.ANNUAL,
                 float_day_count_type: DayCountTypes = DayCountTypes.THIRTY_E_360,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):
        """ Create an overnight index swap contract giving the contract start
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

###############################################################################

    def value(self,
              valuation_date: Date,
              ois_curve: DiscountCurve,
              first_fixing_rate=None):
        """ Value the interest rate swap on a value date given a single Ibor
        discount curve. """

        fixed_leg_value = self._fixed_leg.value(valuation_date,
                                                ois_curve)

        float_leg_value = self._float_leg.value(valuation_date,
                                                ois_curve,
                                                ois_curve,
                                                first_fixing_rate)

        value = fixed_leg_value + float_leg_value
        return value

##########################################################################

    def pv01(self, valuation_date, discount_curve):
        """ Calculate the value of 1 basis point coupon on the fixed leg. """

        pv = self._fixed_leg.value(valuation_date, discount_curve)
        pv01 = pv / self._fixed_leg._coupon / self._fixed_leg._notional

        # Needs to be positive even if it is a payer leg and/or coupon < 0
        pv01 = np.abs(pv01)
        return pv01

##########################################################################

    def swap_rate(self, valuation_date, ois_curve, first_fixing_rate=None):
        """ Calculate the fixed leg coupon that makes the swap worth zero.
        If the valuation date is before the swap payments start then this
        is the forward swap rate as it starts in the future. The swap rate
        is then a forward swap rate and so we use a forward discount
        factor. If the swap fixed leg has begun then we have a spot
        starting swap. """

        pv01 = self.pv01(valuation_date, ois_curve)

        float_leg_value = self._float_leg.value(valuation_date,
                                                ois_curve,
                                                ois_curve,
                                                first_fixing_rate)

        cpn = float_leg_value / pv01 / self._fixed_leg._notional
        return cpn

###############################################################################

    def print_fixed_leg_pv(self):
        """ Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows. """

        self._fixed_leg.print_valuation()

###############################################################################

    def print_float_leg_pv(self):
        """ Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows. """

        self._float_leg.print_valuation()

###############################################################################

    def print_flows(self):
        """ Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows. """

        self._fixed_leg.print_payments()
        self._float_leg.print_payments()

##########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += self._fixed_leg.__repr__()
        s += "\n"
        s += self._float_leg.__repr__()
        return s

###############################################################################

    def _print(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        print(self)

###############################################################################
