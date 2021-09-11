##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.day_count import DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.calendar import CalendarTypes, DateGenRuleTypes
from ...utils.calendar import Calendar, BusDayAdjustTypes
from ...utils.helpers import check_argument_types, label_to_string
from ...utils.math import ONE_MILLION
from ...utils.global_types import SwapTypes
from ...market.curves.discount_curve import DiscountCurve

from .swap_float_leg import SwapFloatLeg

###############################################################################


class IborBasisSwap:
    """ Class for managing an Ibor-Ibor basis swap contract. This is a
    contract in which a floating leg with one LIBOR tenor is exchanged for a
    floating leg payment in a different LIBOR tenor. There is no exchange of
    par. The contract is entered into at zero initial cost. The contract lasts
    from an effective date to a specified maturity date.

    The value of the contract is the NPV of the two coupon streams. Discounting
    is done on a supplied discount curve which can be different from the two
    index discount from which the implied index rates are extracted. """

    def __init__(self,
                 effective_date: Date,  # Date interest starts to accrue
                 termination_date_or_tenor: (Date, str),  # Date contract ends
                 leg1Type: SwapTypes,
                 leg1FreqType: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 leg1DayCountType: DayCountTypes = DayCountTypes.THIRTY_E_360,
                 leg1Spread: float = 0.0,
                 leg2FreqType: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 leg2DayCountType: DayCountTypes = DayCountTypes.THIRTY_E_360,
                 leg2Spread: float = 0.0,
                 notional: float = ONE_MILLION,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):
        """ Create a Ibor basis swap contract giving the contract start
        date, its maturity, frequency and day counts on the two floating 
        legs and notional. The floating leg parameters have default
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

        leg2Type = SwapTypes.PAY
        if leg1Type == SwapTypes.PAY:
            leg2Type = SwapTypes.RECEIVE

        payment_lag = 0
        principal = 0.0

        self._floatLeg1 = SwapFloatLeg(effective_date,
                                       self._termination_date,
                                       leg1Type,
                                       leg1Spread,
                                       leg1FreqType,
                                       leg1DayCountType,
                                       notional,
                                       principal,
                                       payment_lag,
                                       calendar_type,
                                       bus_day_adjust_type,
                                       date_gen_rule_type)

        self._floatLeg2 = SwapFloatLeg(effective_date,
                                       self._termination_date,
                                       leg2Type,
                                       leg2Spread,
                                       leg2FreqType,
                                       leg2DayCountType,
                                       notional,
                                       principal,
                                       payment_lag,
                                       calendar_type,
                                       bus_day_adjust_type,
                                       date_gen_rule_type)

###############################################################################

    def value(self,
              valuation_date: Date,
              discount_curve: DiscountCurve,
              index_curveLeg1: DiscountCurve = None,
              index_curveLeg2: DiscountCurve = None,
              firstFixingRateLeg1=None,
              firstFixingRateLeg2=None):
        """ Value the interest rate swap on a value date given a single Ibor
        discount curve and an index curve for the Ibors on each swap leg. """

        if index_curveLeg1 is None:
            index_curveLeg1 = discount_curve

        if index_curveLeg2 is None:
            index_curveLeg2 = discount_curve

        floatLeg1Value = self._floatLeg1.value(valuation_date,
                                               discount_curve,
                                               index_curveLeg1,
                                               firstFixingRateLeg1)

        floatLeg2Value = self._floatLeg2.value(valuation_date,
                                               discount_curve,
                                               index_curveLeg2,
                                               firstFixingRateLeg2)

        value = floatLeg1Value + floatLeg2Value
        return value

###############################################################################

    def print_float_leg_1_pv(self):
        """ Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows. """

        self._floatLeg1.print_valuation()

###############################################################################

    def print_float_leg_2_pv(self):
        """ Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows. """

        self._floatLeg2.print_valuation()

###############################################################################

    def print_flows(self):
        """ Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows. """

        self._floatLeg1.print_payments()
        self._floatLeg2.print_payments()

##########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += self._floatLeg1.__repr__()
        s += "\n"
        s += self._floatLeg2.__repr__()
        return s

###############################################################################

    def _print(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        print(self)

###############################################################################
