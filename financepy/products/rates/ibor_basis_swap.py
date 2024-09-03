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
    """Class for managing an Ibor-Ibor basis swap contract. This is a
    contract in which a floating leg with one LIBOR tenor is exchanged for a
    floating leg payment in a different LIBOR tenor. There is no exchange of
    par. The contract is entered into at zero initial cost. The contract lasts
    from an effective date to a specified maturity date.

    The value of the contract is the NPV of the two coupon streams. Discounting
    is done on a supplied discount curve which can be different from the two
    index discount from which the implied index rates are extracted."""

    def __init__(
        self,
        effective_dt: Date,  # Date interest starts to accrue
        term_dt_or_tenor: (Date, str),  # Date contract ends
        leg_1_type: SwapTypes,
        leg_1_freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
        leg_1_day_count_type: DayCountTypes = DayCountTypes.THIRTY_E_360,
        leg_1_spread: float = 0.0,
        leg2FreqType: FrequencyTypes = FrequencyTypes.QUARTERLY,
        leg2DayCountType: DayCountTypes = DayCountTypes.THIRTY_E_360,
        leg2Spread: float = 0.0,
        notional: float = ONE_MILLION,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
    ):
        """Create a Ibor basis swap contract giving the contract start
        date, its maturity, frequency and day counts on the two floating
        legs and notional. The floating leg parameters have default
        values that can be overwritten if needed. The start date is contractual
        and is the same as the settlement date for a new swap. It is the date
        on which interest starts to accrue. The end of the contract is the
        termination date. This is not adjusted for business days. The adjusted
        termination date is called the maturity date. This is calculated."""

        check_argument_types(self.__init__, locals())

        if isinstance(term_dt_or_tenor, Date):
            self.termination_dt = term_dt_or_tenor
        else:
            self.termination_dt = effective_dt.add_tenor(term_dt_or_tenor)

        calendar = Calendar(cal_type)
        self.maturity_dt = calendar.adjust(self.termination_dt, bd_type)

        if effective_dt > self.maturity_dt:
            raise FinError("Start date after maturity date")

        leg2Type = SwapTypes.PAY
        if leg_1_type == SwapTypes.PAY:
            leg2Type = SwapTypes.RECEIVE

        payment_lag = 0
        principal = 0.0

        self.float_leg_1 = SwapFloatLeg(
            effective_dt,
            self.termination_dt,
            leg_1_type,
            leg_1_spread,
            leg_1_freq_type,
            leg_1_day_count_type,
            notional,
            principal,
            payment_lag,
            cal_type,
            bd_type,
            dg_type,
        )

        self.float_leg_2 = SwapFloatLeg(
            effective_dt,
            self.termination_dt,
            leg2Type,
            leg2Spread,
            leg2FreqType,
            leg2DayCountType,
            notional,
            principal,
            payment_lag,
            cal_type,
            bd_type,
            dg_type,
        )

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        discount_curve: DiscountCurve,
        index_curve_leg_1: DiscountCurve = None,
        index_curve_leg_2: DiscountCurve = None,
        first_fixing_rate_leg_1=None,
        first_fixing_rate_leg_2=None,
    ):
        """Value the interest rate swap on a value date given a single Ibor
        discount curve and an index curve for the Ibors on each swap leg."""

        if index_curve_leg_1 is None:
            index_curve_leg_1 = discount_curve

        if index_curve_leg_2 is None:
            index_curve_leg_2 = discount_curve

        float_leg_1Value = self.float_leg_1.value(
            value_dt,
            discount_curve,
            index_curve_leg_1,
            first_fixing_rate_leg_1,
        )

        float_leg_2Value = self.float_leg_2.value(
            value_dt,
            discount_curve,
            index_curve_leg_2,
            first_fixing_rate_leg_2,
        )

        value = float_leg_1Value + float_leg_2Value
        return value

    ###########################################################################

    def print_float_leg_1_pv(self):
        """Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows."""

        self.float_leg_1.print_valuation()

    ###########################################################################

    def print_float_leg_2_pv(self):
        """Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows."""

        self.float_leg_2.print_valuation()

    ###########################################################################

    def print_payments(self):
        """Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows."""

        self.float_leg_1.print_payments()
        self.float_leg_2.print_payments()

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += self.float_leg_1.__repr__()
        s += "\n"
        s += self.float_leg_2.__repr__()
        return s

    ###########################################################################

    def _print(self):
        """Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond."""
        print(self)


###############################################################################
