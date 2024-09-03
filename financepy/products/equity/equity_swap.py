##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.day_count import DayCountTypes
from ...utils.frequency import FrequencyTypes, annual_frequency
from ...utils.calendar import CalendarTypes, DateGenRuleTypes
from ...utils.calendar import Calendar, BusDayAdjustTypes
from ...utils.helpers import check_argument_types
from ...utils.global_types import SwapTypes, ReturnTypes
from ...market.curves.discount_curve import DiscountCurve
from ...products.rates.swap_float_leg import SwapFloatLeg
from ...products.equity.equity_swap_leg import EquitySwapLeg

###############################################################################


class EquitySwap:
    """Class for managing a standard Equity vs Float leg swap. This is a
    contract in which an equity payment leg is exchanged for a series of
    floating rates payments. There is no exchange of principal. The contract
    is entered into at zero initial cost when spreads are zero. The contract
    lasts from an effective date to a specified maturity date.

    The equity payments are not known fully until the end of payment period.

    The floating rate is not known fully until the end of the preceding payment
    period. It is set in advance and paid in arrears.

    The value of the contract is the NPV of the two coupon streams. Discounting
    is done on a supplied discount curve which is separate from the curve from
    which the implied index rates are extracted."""

    def __init__(
        self,
        effective_dt: Date,  # Date contract starts or last Equity Reset
        term_dt_or_tenor: (Date, str),  # Date contract ends
        eq_leg_type: SwapTypes,
        eq_freq_type: FrequencyTypes,
        eq_dc_type: DayCountTypes,
        strike: float,  # Price at effective date
        quantity: float = 1.0,  # Quantity at effective date
        eq_payment_lag: int = 0,
        eq_return_type: ReturnTypes = ReturnTypes.TOTAL_RETURN,
        rate_freq_type: FrequencyTypes = FrequencyTypes.MONTHLY,
        rate_dc_type: DayCountTypes = DayCountTypes.ACT_360,
        rate_spread: float = 0.0,
        rate_payment_lag: int = 0,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
        end_of_month: bool = False,
    ):
        """Create an equity swap contract given the contract effective date, its maturity,
        underlying price and quantity, day count convention and return type and other details.
        The equity leg parameters have default values that can be overwritten if needed.
        The start date is contractual and is the same as the settlement date for a new swap.
        It is the date on which interest starts to accrue. The end of the contract is the
        termination date. This is not adjusted for business days. The adjusted termination
        date is called the maturity date. This is calculated."""

        check_argument_types(self.__init__, locals())

        if isinstance(term_dt_or_tenor, Date):
            self.termination_dt = term_dt_or_tenor
        else:
            self.termination_dt = effective_dt.add_tenor(term_dt_or_tenor)

        calendar = Calendar(cal_type)
        self.maturity_dt = calendar.adjust(self.termination_dt, bd_type)

        if effective_dt > self.maturity_dt:
            raise FinError("Start date after maturity date")

        self.effective_dt = effective_dt

        # There is no exchange of principal
        self.principal = 0.0

        rate_leg_type = SwapTypes.PAY
        if eq_leg_type == SwapTypes.PAY:
            rate_leg_type = SwapTypes.RECEIVE

        self.equity_leg = EquitySwapLeg(
            effective_dt,
            self.maturity_dt,
            eq_leg_type,
            eq_freq_type,
            eq_dc_type,
            strike,
            quantity,
            eq_payment_lag,
            eq_return_type,
            cal_type,
            bd_type,
            dg_type,
            end_of_month,
        )

        # Fixed Rate Leg not implemented yet
        self.rate_leg = SwapFloatLeg(
            effective_dt,
            self.maturity_dt,
            rate_leg_type,
            rate_spread,
            rate_freq_type,
            rate_dc_type,
            self.equity_leg.notional,
            self.principal,
            rate_payment_lag,
            cal_type,
            bd_type,
            dg_type,
            end_of_month,
        )

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        discount_curve: DiscountCurve,
        index_curve: DiscountCurve = None,
        dividend_curve: DiscountCurve = None,
        current_price: float = None,
        first_fixing_rate=None,
    ):
        """Value the Equity swap on a valuation date."""

        self.equity_leg_value = self.equity_leg.value(
            value_dt,
            discount_curve,
            index_curve,
            dividend_curve,
            current_price,
        )
        self._fill_rate_notional_array()

        self.rate_leg_value = self.rate_leg.value(
            value_dt, discount_curve, index_curve, first_fixing_rate
        )

        return self.equity_leg_value + self.rate_leg_value

    ###########################################################################

    def _fill_rate_notional_array(self):
        """In an equity swap, at every equity reset, the notional
        of the contract is updated to reflect the new underlying
        price.

        This is a helper function that takes the Equity Notional list
        from Equity Leg and convert it to a Notional array that fits
        the payment schedule defined for the rate leg.
        """

        # Assumption: Rate frequency type is a multiple of Equity's
        eq_freq = annual_frequency(self.equity_leg.freq_type)
        rate_freq = annual_frequency(self.rate_leg.freq_type)

        multiple = int(rate_freq // eq_freq)
        is_multiple = int(rate_freq % eq_freq) == 0

        if (eq_freq is None or rate_freq is None) or (not is_multiple):
            raise FinError("Invalid frequency type assigned!")

        self.rate_leg.notional_array = []
        for last_notional in self.equity_leg.last_notionals:
            for _ in range(multiple):
                self.rate_leg.notional_array.append(last_notional)

    ###########################################################################

    def __repr__(self):
        s = "EQUITY LEG:\n"
        s += self.equity_leg.__repr__()
        s += "\nRATE LEG:\n"
        s += self.rate_leg.__repr__()

        return s

    ###########################################################################
