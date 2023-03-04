##############################################################################
# Copyright (C) 2023 Dominic O'Kane
##############################################################################

from enum import Enum

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.day_count import DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.calendar import CalendarTypes, DateGenRuleTypes
from ...utils.calendar import Calendar, BusDayAdjustTypes
from ...utils.helpers import check_argument_types
from ...utils.global_types import SwapTypes, ReturnTypes
from ...market.curves.discount_curve import DiscountCurve
from ...products.rates.swap_fixed_leg import SwapFixedLeg
from ...products.rates.swap_float_leg import SwapFloatLeg
from ...products.equity.equity_swap_leg import SwapEquityLeg


###############################################################################
class EquitySwap:
    """ Class for managing a standard Equity vs Float leg swap. This is a contract 
    in which an equity payment leg is exchanged for a series of either floating rates 
    payments. There is no exchange of principal. The contract is entered into at zero 
    initial cost. The contract lasts from a start date to a specified maturity date.

    The equity payments  are not known fully until the end of the valuation period. 

    The value of the contract is the NPV of the two coupon streams. Discounting
    is done on a supplied discount curve which is separate from the curve from
    which the implied index rates are extracted. """

    def __init__(self,
                 effective_date: Date, ## Date contract starts or last Equity Reset
                 termination_date_or_tenor: (Date, str), ## Date contract ends
                 eq_leg_type: SwapTypes,
                 eq_freq_type: FrequencyTypes,
                 eq_day_count_type: DayCountTypes,
                 strike: float, ## Price at effective date
                 quantity: float = 1.0, ## Quantity at effective date
                 eq_payment_lag: int = 0,
                 eq_return_type: ReturnTypes = ReturnTypes.TOTAL_RETURN,
                 rate_freq_type: FrequencyTypes = FrequencyTypes.ANNUAL, 
                 rate_day_count_type: DayCountTypes = DayCountTypes.ACT_360,
                 rate_spread: float = 0.0, 
                 rate_payment_lag: int = 0,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
                 end_of_month: bool = False):
        """ Create an equity swap contract giving the contract effective
        date, its maturity, notional based on underlying price & quantity, 
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

        ## There is no exchange of principal
        self._principal = 0.0

        rate_leg_type = SwapTypes.PAY
        if eq_leg_type == SwapTypes.PAY:
            rate_leg_type = SwapTypes.RECEIVE

        self._equity_leg = SwapEquityLeg(effective_date,
                                            self._maturity_date,
                                            eq_leg_type,
                                            eq_freq_type,
                                            eq_day_count_type,
                                            strike,
                                            quantity,
                                            eq_payment_lag,
                                            eq_return_type,
                                            calendar_type,
                                            bus_day_adjust_type,
                                            date_gen_rule_type,
                                            end_of_month)

        ## Fixed Rate Leg not implemented yet
        self._rate_leg = SwapFloatLeg(effective_date,
                                            self._maturity_date,
                                            rate_leg_type,
                                            rate_spread,
                                            rate_freq_type,
                                            rate_day_count_type,
                                            self._equity_leg._notional,
                                            self._principal,
                                            rate_payment_lag,
                                            calendar_type,
                                            bus_day_adjust_type,
                                            date_gen_rule_type,
                                            end_of_month)

    ###########################################################################

    def value(self,
              valuation_date: Date,
              discount_curve: DiscountCurve,
              index_curve: DiscountCurve = None,
              dividend_curve: DiscountCurve = None,
              current_price: float = None,
              firstFixingRate=None):
        """ Value the Equity swap on a valuation date. """

        self._equity_leg_value = self._equity_leg.value(valuation_date,
                                                        discount_curve,
                                                        index_curve,
                                                        dividend_curve,
                                                        current_price,
                                                        firstFixingRate)
        
        ## Notionals from equity leg are computed and can fill rate leg   
        self._rate_leg._notional_array = self._equity_leg._last_notionals

        self._rate_leg_value = self._rate_leg.value(valuation_date,
                                                    discount_curve,
                                                    index_curve,
                                                    firstFixingRate)
        
        return self._equity_leg_value + self._rate_leg_value