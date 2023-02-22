##############################################################################
# Copyright (C) 2023 Dominic O'Kane
##############################################################################


from enum import Enum

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.global_vars import gSmall
from ...utils.day_count import DayCountTypes, 
from ...utils.frequency import FrequencyTypes, annual_frequency
from ...utils.calendar import CalendarTypes, DateGenRuleTypes
from ...utils.calendar import Calendar, BusDayAdjustTypes
from ...utils.helpers import check_argument_types, label_to_string
from ...utils.global_types import SwapTypes, ReturnTypes, RateTypes
from ...market.curves.discount_curve import DiscountCurve
from ...products.rates.swap_fixed_leg import SwapFixedLeg
from ...products.rates.swap_float_leg import SwapFloatLeg

from equity_swap_leg import SwapEquityLeg


###############################################################################
class EquitySwap:
    """ Class for managing a standard Equity vs fixed-or-float leg swap. This is 
    a contract in which a equity payment leg is exchanged for a series of either
    fixed or floating rates payments. There is no exchange of principal. The 
    contract is entered into at zero initial cost. The contract lasts from
    a start date to a specified maturity date.

    Neither the equity payments nor the floating rate are not known fully until 
    the end of the preceding payment period. It is set in advance and paid in arrears.

    The value of the contract is the NPV of the two coupon streams. Discounting
    is done on a supplied discount curve which is separate from the curve from
    which the implied index rates are extracted. """

    def __init__(self,
                 effective_date: Date,  # Date equity effective starts
                 maturity_date_or_tenor: (Date, str),  #  Date contract ends

                 ## Equity args
                 leg_type: SwapTypes,
                 eq_freq_type: FrequencyTypes,
                 eq_day_count_type: DayCountTypes,
                 underlyer_price: float,  # Price at effective date
                 underlyer_quantity: int = 1, # Quantity at effective date
                 eq_payment_lag: int = 0,
                 eq_return_type: ReturnTypes = ReturnTypes.TOTAL_RETURN,

                 ### Rate args
                 rate_leg_type: RateTypes = RateTypes.FIXED,
                 rate_freq_type: FrequencyTypes = FrequencyTypes.ANNUAL, 
                 rate_day_count_type: DayCountTypes = DayCountTypes.ACT_360,
                 rate_spread: float = 0.0, ## Applied if float
                 rate_payment_lag: int = 0,

                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
                 end_of_month: bool = False):
        """ Create an equity swap contract giving the contract start
        date, its maturity, fixed coupon, fixed leg frequency, fixed leg day
        count convention and notional. The floating leg parameters have default
        values that can be overwritten if needed. The start date is contractual
        and is the same as the settlement date for a new swap. It is the date
        on which interest starts to accrue. The end of the contract is the
        termination date. This is not adjusted for business days. The adjusted
        termination date is called the maturity date. This is calculated. """

        check_argument_types(self.__init__, locals())

        if type(maturity_date_or_tenor) == Date:
            self._maturity_date = maturity_date_or_tenor
        else:
            self._maturity_date = effective_date.add_tenor(
                maturity_date_or_tenor)

        calendar = Calendar(calendar_type)
        self._maturity_date = calendar.adjust(self._maturity_date,
                                              bus_day_adjust_type)

        if effective_date > self._maturity_date:
            raise FinError("Start date after maturity date")

        self._effective_date = effective_date

        self._equity_leg = SwapEquityLeg(effective_date,
        )

        if rate_leg_type == RateTypes.FIXED:
            self._rate_leg = SwapFixedLeg()
        else:
            self._rate_leg = SwapFloatLeg()

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

    