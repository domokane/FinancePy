##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
import pandas as pd

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.global_vars import g_small
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
    """Class for managing a standard Fixed vs IBOR swap. This is a contract
    in which a fixed payment leg is exchanged for a series of floating rates
    payments linked to some IBOR index rate. There is no exchange of principal.
    The contract is entered into at zero initial cost. The contract lasts from
    a start date to a specified maturity date.

    The floating rate is not known fully until the end of the preceding payment
    period. It is set in advance and paid in arrears.

    The value of the contract is the NPV of the two cpn streams. Discounting
    is done on a supplied discount curve which is separate from the curve from
    which the implied index rates are extracted."""

    def __init__(
        self,
        effective_dt: Date,  # Date interest starts to accrue
        term_dt_or_tenor: (Date, str),  # Date contract ends
        fixed_leg_type: SwapTypes,
        fixed_cpn: float,  # Fixed cpn (annualised)
        fixed_freq_type: FrequencyTypes,
        fixed_dc_type: DayCountTypes,
        notional: float = ONE_MILLION,
        float_spread: float = 0.0,
        float_freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
        float_dc_type: DayCountTypes = DayCountTypes.THIRTY_E_360,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
    ):
        """Create an interest rate swap contract giving the contract start
        date, its maturity, fixed cpn, fixed leg frequency, fixed leg day
        count convention and notional. The floating leg parameters have default
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

        self.effective_dt = effective_dt

        float_leg_type = SwapTypes.PAY
        if fixed_leg_type == SwapTypes.PAY:
            float_leg_type = SwapTypes.RECEIVE

        payment_lag = 0
        principal = 0.0

        self.fixed_leg = SwapFixedLeg(
            effective_dt,
            self.termination_dt,
            fixed_leg_type,
            fixed_cpn,
            fixed_freq_type,
            fixed_dc_type,
            notional,
            principal,
            payment_lag,
            cal_type,
            bd_type,
            dg_type,
        )

        self.float_leg = SwapFloatLeg(
            effective_dt,
            self.termination_dt,
            float_leg_type,
            float_spread,
            float_freq_type,
            float_dc_type,
            notional,
            principal,
            payment_lag,
            cal_type,
            bd_type,
            dg_type,
        )

    ###########################################################################

    def get_fixed_rate(self):
        """
        easy read access to the coupon (fixed rate)
        """
        return self.fixed_leg.cpn

    ###########################################################################

    def set_fixed_rate(self, new_rate: float):
        """
        Sometimes we need to reset the coupon (fixed rate)
        This function updates caches that depend on it
        """
        self.fixed_leg.cpn = new_rate
        self.fixed_leg.generate_payments()

    ###########################################################################

    def set_fixed_rate_to_atm(
        self,
        valuation_date: Date,
        discount_curve: DiscountCurve,
        index_curve: DiscountCurve = None,
        first_fixing: float = None,
    ):
        """
        Reset fixed rate to atm given curve(s). returns the new atm
        """
        atm = self.swap_rate(
            valuation_date, discount_curve, index_curve, first_fixing
        )
        self.set_fixed_rate(atm)
        return atm

    ###########################################################################
    def value(
        self,
        value_dt: Date,
        discount_curve: DiscountCurve,
        index_curve: DiscountCurve = None,
        first_fixing_rate=None,
        pv_only=True,
    ):
        """Value the interest rate swap on a value date given a single Ibor
        discount curve."""

        if index_curve is None:
            index_curve = discount_curve

        fixed_leg_results = self.fixed_leg.value(
            value_dt, discount_curve, pv_only=pv_only
        )

        float_leg_results = self.float_leg.value(
            value_dt,
            discount_curve,
            index_curve,
            first_fixing_rate,
            pv_only=pv_only,
        )

        if pv_only:
            value = fixed_leg_results + float_leg_results
            return value
        else:
            value = fixed_leg_results[0] + float_leg_results[0]
            cashflow_report = pd.concat(
                [fixed_leg_results[1], float_leg_results[1]], ignore_index=True
            )

            return value, cashflow_report

    ###########################################################################

    def valuation_details(
        self,
        valuation_date: Date,
        discount_curve: DiscountCurve,
        index_curve: DiscountCurve = None,
        firstFixingRate=None,
    ):
        """
        A long-hand method that returns various details relevant to valuation in a dictionary
        Slower than value(...) so should not be used when performance is important

        We want the output dictionary to have  the same labels for different bechmarks
        (depos, fras, swaps) because we want to present them together so please do not stick new outputs into
        one of them only
        """
        if index_curve is None:
            index_curve = discount_curve

        fixed_leg_value = self.fixed_leg.value(valuation_date, discount_curve)

        float_leg_value = self.float_leg.value(
            valuation_date, discount_curve, index_curve, firstFixingRate
        )

        value = fixed_leg_value + float_leg_value
        pv01 = np.abs(
            fixed_leg_value / self.fixed_leg.cpn / self.fixed_leg.notional
        )
        pay_receive_float = (
            -1 if self.float_leg.leg_type == SwapTypes.PAY else 1
        )
        swap_rate = (
            float_leg_value
            / self.float_leg.notional
            / pv01
            / pay_receive_float
        )

        # VP: There is significant amount of confusion here with swap_type vs notional.
        is_payers = (
            self.fixed_leg.leg_type == SwapTypes.PAY
            and self.fixed_leg.notional > 0
        ) or (
            self.fixed_leg.leg_type == SwapTypes.RECEIVE
            and self.fixed_leg.notional < 0
        )

        pvbp_sign = 1 if is_payers else -1

        out = {
            "type": type(self).__name__,
            "start_date": self.effective_dt,
            "maturity_date": self.maturity_dt,
            "day_count_type": self.fixed_leg.dc_type.name,
            "fixed_leg_type": self.fixed_leg.leg_type.name,
            "fixed_freq_type": self.fixed_leg.freq_type.name,
            "notional": self.fixed_leg.notional,
            "contract_rate": self.fixed_leg.cpn,
            "market_rate": swap_rate,
            "spot_pvbp": pv01 * pvbp_sign,
            "fwd_pvbp": pv01
            * pvbp_sign
            / discount_curve.df(self.effective_dt),
            "unit_value": value / self.fixed_leg.notional,
            "value": value,
            # ignoring bus day adj type, calendar, etc for now
        }
        return out

    ###########################################################################

    def pv01(self, value_dt, discount_curve):
        """Calculate the value of 1 basis point coupon on the fixed leg."""

        pv = self.fixed_leg.value(value_dt, discount_curve)
        pv01 = pv / self.fixed_leg.cpn / self.fixed_leg.notional
        # Needs to be positive even if it is a payer leg
        pv01 = np.abs(pv01)
        return pv01

    ###########################################################################

    def swap_rate(
        self,
        value_dt: Date,
        discount_curve: DiscountCurve,
        index_curve: DiscountCurve = None,
        first_fixing: float = None,
    ):
        """Calculate the fixed leg cpn that makes the swap worth zero.
        If the valuation date is before the swap payments start then this
        is the forward swap rate as it starts in the future. The swap rate
        is then a forward swap rate and so we use a forward discount
        factor. If the swap fixed leg has begun then we have a spot
        starting swap. The swap rate can also be calculated in a dual curve
        approach but in this case the first fixing on the floating leg is
        needed."""

        pv01 = self.pv01(value_dt, discount_curve)

        if abs(pv01) < g_small:
            raise FinError("PV01 is zero. Cannot compute swap rate.")

        # VP: I commented out this shortcut below because it is inconsistent with value(...) function
        # there are some subtle differences due to day_counts
        """
        float_leg_pv = 0.0

        if valuation_date < self.effective_dt:
            df0 = discount_curve.df(self.effective_dt)
        else:
            df0 = discount_curve.df(value_dt)

        if index_curve is None:
            df_t = discount_curve.df(self.maturity_dt)
            float_leg_pv = (df0 - df_t)

        else:
            float_leg_pv = self.float_leg.value(value_dt,
                                                discount_curve,
                                                index_curve,
                                                first_fixing)

            float_leg_pv /= self.fixed_leg.notional
        """
        # VP: this is more consistent with value(..):
        float_leg_pv = self.float_leg.value(
            value_dt, discount_curve, index_curve, first_fixing
        )

        float_leg_pv /= self.float_leg.notional

        # Make sure we get the sign right
        if self.float_leg.leg_type == SwapTypes.PAY:
            float_leg_pv = -float_leg_pv

        cpn = float_leg_pv / pv01
        return cpn

    ###########################################################################

    def cash_settled_pv01(self, value_dt, flat_swap_rate, freq_type):
        """Calculate the forward value of an annuity of a forward starting
        swap using a single flat discount rate equal to the swap rate. This is
        used in the pricing of a cash-settled swaption in the IborSwaption
        class. This method does not affect the standard valuation methods."""

        m = annual_frequency(freq_type)

        if m == 0:
            raise FinError("Frequency cannot be zero.")

        """ The swap may have started in the past but we can only value
        payments that have occurred after the valuation date. """
        start_index = 0
        while self.fixed_leg.payment_dts[start_index] < value_dt:
            start_index += 1

        """ If the swap has yet to settle then we do not include the
        start date of the swap as a cpn payment date. """
        if value_dt <= self.effective_dt:
            start_index = 1

        """ Now PV fixed leg flows. """
        flat_pv01 = 0.0
        df = 1.0
        alpha = 1.0 / m

        for _ in self.fixed_leg.payment_dts[start_index:]:
            df = df / (1.0 + alpha * flat_swap_rate)
            flat_pv01 += df * alpha

        return flat_pv01

    ###########################################################################

    def print_fixed_leg_pv(self):
        """Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows."""

        self.fixed_leg.print_valuation()

    ###########################################################################

    def print_float_leg_pv(self):
        """Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows."""

        self.float_leg.print_valuation()

    ###########################################################################

    def print_payments(self):
        """Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows."""

        self.fixed_leg.print_payments()
        self.float_leg.print_payments()

    ###########################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += self.fixed_leg.__repr__()
        s += "\n"
        s += self.float_leg.__repr__()
        return s

    ###########################################################################

    def _print(self):
        """Print a list of the unadjusted cpn payment dates used in
        analytic calculations for the bond."""
        print(self)


###############################################################################
