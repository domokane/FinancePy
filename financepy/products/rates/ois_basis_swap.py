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


class OISBasisSwap:
    """ Class for managing an Ibor-OIS basis swap contract. This is a
    contract in which a floating leg with one LIBOR tenor is exchanged for a
    floating leg payment of an overnight index swap. There is no exchange of
    par. The contract is entered into at zero initial cost. The contract lasts
    from a start date to a specified maturity date.

    The value of the contract is the NPV of the two coupon streams. Discounting
    is done on a supplied discount curve which is separate from the discount
    from which the implied index rates are extracted. """

    def __init__(self,
                 effective_dt: Date,  # Date interest starts to accrue
                 term_dt_or_tenor: (Date, str),  # Date contract ends
                 ibor_type: SwapTypes,
                 ibor_freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 ibor_day_count_type: DayCountTypes = DayCountTypes.THIRTY_E_360,
                 ibor_spread: float = 0.0,
                 ois_freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 ois_day_count_type: DayCountTypes = DayCountTypes.THIRTY_E_360,
                 ois_spread: float = 0.0,
                 ois_payment_lag: int = 0,
                 notional: float = ONE_MILLION,
                 cal_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):
        """ Create a Ibor basis swap contract giving the contract start
        date, its maturity, frequency and day counts on the two floating
        legs and notional. The floating leg parameters have default
        values that can be overwritten if needed. The start date is contractual
        and is the same as the settlement date for a new swap. It is the date
        on which interest starts to accrue. The end of the contract is the
        termination date. This is not adjusted for business days. The adjusted
        termination date is called the maturity date. This is calculated. """

        check_argument_types(self.__init__, locals())

        if isinstance(term_dt_or_tenor, Date):
            self.termination_dt = term_dt_or_tenor
        else:
            self.termination_dt = effective_dt.add_tenor(
                term_dt_or_tenor)

        calendar = Calendar(cal_type)
        self.maturity_dt = calendar.adjust(self.termination_dt,
                                           bd_type)

        if effective_dt > self.maturity_dt:
            raise FinError("Start date after maturity date")

        ois_type = SwapTypes.PAY
        if ibor_type == SwapTypes.PAY:
            ois_type = SwapTypes.RECEIVE

        principal = 0.0

        self.float_ibor_leg = SwapFloatLeg(effective_dt,
                                           self.termination_dt,
                                           ibor_type,
                                           ibor_spread,
                                           ibor_freq_type,
                                           ibor_day_count_type,
                                           notional,
                                           principal,
                                           0,
                                           cal_type,
                                           bd_type,
                                           dg_type)

        self.float_ois_leg = SwapFloatLeg(effective_dt,
                                          self.termination_dt,
                                          ois_type,
                                          ois_spread,
                                          ois_freq_type,
                                          ois_day_count_type,
                                          notional,
                                          principal,
                                          ois_payment_lag,
                                          cal_type,
                                          bd_type,
                                          dg_type)

###############################################################################

    def value(self,
              value_dt: Date,
              discount_curve: DiscountCurve,
              index_ibor_curve: DiscountCurve = None,
              index_ois_curve: DiscountCurve = None,
              first_fixing_rate_leg_1=None,
              first_fixing_rate_leg_2=None):
        """ Value the interest rate swap on a value date given a single Ibor
        discount curve and an index curve for the Ibors on each swap leg. """

        if index_ibor_curve is None:
            index_ibor_curve = discount_curve

        if index_ois_curve is None:
            index_ois_curve = discount_curve

        float_ibor_leg_value = self.float_ibor_leg.value(value_dt,
                                                         discount_curve,
                                                         index_ibor_curve,
                                                         first_fixing_rate_leg_1)

        float_ois_leg_value = self.float_ois_leg.value(value_dt,
                                                       discount_curve,
                                                       index_ois_curve,
                                                       first_fixing_rate_leg_2)

        value = float_ibor_leg_value + float_ois_leg_value
        return value

###############################################################################

    def print_payments(self):
        """ Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows. """

        self.float_ibor_leg.print_payments()
        self.float_ois_leg.print_payments()

##########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += self.float_ibor_leg._repr__()
        s += "\n"
        s += self.float_ois_leg._repr__()
        return s

###############################################################################

    def _print(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        print(self)

###############################################################################
