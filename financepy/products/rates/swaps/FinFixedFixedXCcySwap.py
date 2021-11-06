##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...utils.FinError import FinError
from ...utils.Date import Date
from ...utils.FinGlobalVariables import gSmall
from ...utils.FinDayCount import FinDayCount, DayCountTypes
from ...utils.FinFrequency import FrequencyTypes, FinFrequency
from ...utils.FinCalendar import CalendarTypes,  DateGenRuleTypes
from ...utils.FinCalendar import FinCalendar, BusDayAdjustTypes
from ...utils.FinSchedule import FinSchedule
from ...utils.FinHelperFunctions import label_to_string, check_argument_types
from ...utils.FinMath import ONE_MILLION
from ...utils.FinGlobalTypes import SwapTypes

##########################################################################


class FinFixedFixedXCcySwap
  """ Class for managing a cross currency swap contract. This is a contract
    in which a fixed or floating payment leg in one currency is exchanged for a
    series of fixed or floating rates in a second currency. There is an
    exchange of par. The contract is entered into at zero initial cost and it
    lasts from a start date to a specified maturity date.

    The value of the contract is the NPV of the two coupon streams. Discounting
    is done on a supplied discount discount (one for each leg) which is separate
    from the curve from which the implied index rates are extracted. """

   def __init__(self,
                 effective_date: Date,  # Date interest starts to accrue
                 termination_date_or_tenor: (Date, str),  # Date contract ends
                 fixed_leg_type: SwapTypes,
                 fixed_coupon: float,  # Fixed coupon (annualised)
                 fixed_freq_type: FrequencyTypes,
                 fixed_day_count_type: DayCountTypes,
                 float_spread: float = 0.0,
                 float_freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 float_day_count_type: DayCountTypes = DayCountTypes.THIRTY_E_360,
                 notional: float = ONE_MILLION,
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
        self._notional = notional

        self._fixed_coupon = fixed_coupon
        self._float_spread = float_spread

        self._fixed_frequency_type = fixed_freq_type
        self._float_frequency_type = float_freq_type

        self._fixed_day_count_type = fixed_day_count_type
        self._float_day_count_type = float_day_count_type

        self._fixed_leg_type = fixed_leg_type

        self._calendar_type = calendar_type
        self._bus_day_adjust_type = bus_day_adjust_type
        self._date_gen_rule_type = date_gen_rule_type

        # These are generated immediately as they are for the entire
        # life of the swap. Given a valuation date we can determine
        # which cash flows are in the future and value the swap
        self._generate_fixed_leg_payment_dates()
        self._generate_float_leg_payment_dates()

        self._adjustedMaturityDate = self._adjustedFixedDates[-1]

        # Need to know latest payment date for bootstrap - DO I NEED THIS ??!
        self._lastPaymentDate = self._maturity_date
        if self._adjustedFixedDates[-1] > self._lastPaymentDate:
            self._lastPaymentDate = self._adjustedFixedDates[-1]

        if self._adjustedFloatDates[-1] > self._lastPaymentDate:
            self._lastPaymentDate = self._adjustedFloatDates[-1]

        # NOT TO BE PRINTED
        self._float_year_fracs = []
        self._floatFlows = []
        self._floatFlowPVs = []
        self._floatDfs = []

        self._fixed_year_fracs = []
        self._fixedFlows = []
        self._fixedDfs = []
        self._fixedFlowPVs = []

        self._firstFixingRate = None
        self._valuation_date = None
        self._fixedStartIndex = None

        self._calc_fixed_leg_flows()

##########################################################################

    def value(self,
              valuation_date,
              discount_curve,
              index_curve,
              first_fixing_rate=None,
              principal=0.0):
        """ Value the interest rate swap on a value date given a single Ibor
        discount curve. """

        fixed_leg_value = self.fixed_leg_value(valuation_date,
                                               discount_curve,
                                               principal)

        float_leg_value = self.float_leg_value(valuation_date,
                                               discount_curve,
                                               index_curve,
                                               first_fixing_rate,
                                               principal)

        value = fixed_leg_value - float_leg_value

        if self._fixed_leg_type == SwapTypes.PAY:
            value = value * (-1.0)

        return value

##########################################################################

    def _generate_fixed_leg_payment_dates(self):
        """ Generate the fixed leg payment dates all the way back to
        the start date of the swap which may precede the valuation date"""
        self._adjustedFixedDates = FinSchedule(
            self._effective_date,
            self._termination_date,
            self._fixed_frequency_type,
            self._calendar_type,
            self._bus_day_adjust_type,
            self._date_gen_rule_type)._generate()

##########################################################################

    def _generate_float_leg_payment_dates(self):
        """ Generate the floating leg payment dates all the way back to
        the start date of the swap which may precede the valuation date"""
        self._adjustedFloatDates = FinSchedule(
            self._effective_date,
            self._termination_date,
            self._float_frequency_type,
            self._calendar_type,
            self._bus_day_adjust_type,
            self._date_gen_rule_type)._generate()

##########################################################################

    def fixed_dates(self):
        """ return a vector of the fixed leg payment dates """
        if self._adjustedFixedDates is None:
            raise FinError("Fixed dates have not been generated")

        return self._adjustedFixedDates[1:]

##########################################################################

    def float_dates(self):
        """ return a vector of the fixed leg payment dates """
        if self._adjustedFloatDates is None:
            raise FinError("Float dates have not been generated")

        return self._adjustedFloatDates[1:]

##########################################################################

    def pv01(self, valuation_date, discount_curve):
        """ Calculate the value of 1 basis point coupon on the fixed leg. """

        pv = self.fixed_leg_value(valuation_date, discount_curve)
        pv01 = pv / self._fixed_coupon / self._notional
        return pv01

##########################################################################

    def swap_rate(self, valuation_date, discount_curve):
        """ Calculate the fixed leg coupon that makes the swap worth zero.
        If the valuation date is before the swap payments start then this
        is the forward swap rate as it starts in the future. The swap rate
        is then a forward swap rate and so we use a forward discount
        factor. If the swap fixed leg has begun then we have a spot
        starting swap. """

        pv01 = self.pv01(valuation_date, discount_curve)

        if valuation_date < self._effective_date:
            df0 = discount_curve.df(self._effective_date)
        else:
            df0 = discount_curve.df(valuation_date)

        dfT = discount_curve.df(self._maturity_date)

        if abs(pv01) < gSmall:
            raise FinError("PV01 is zero. Cannot compute swap rate.")

        cpn = (df0 - dfT) / pv01
        return cpn

##########################################################################

    def fixed_leg_value(self, valuation_date, discount_curve, principal=0.0):

        self._valuation_date = valuation_date

        self._fixed_year_fracs = []
        self._fixedFlows = []
        self._fixedDfs = []
        self._fixedFlowPVs = []
        self._fixed_total_pv = []

        day_counter = FinDayCount(self._fixed_day_count_type)

        """ The swap may have started in the past but we can only value
        payments that have occurred after the valuation date. """
        start_index = 0
        while self._adjustedFixedDates[start_index] < valuation_date:
            start_index += 1

        """ If the swap has yet to settle then we do not include the
        start date of the swap as a coupon payment date. """
        if valuation_date <= self._effective_date:
            start_index = 1

        self._fixedStartIndex = start_index

        """ Now PV fixed leg flows. """
        self._dfValuationDate = discount_curve.df(valuation_date)

        pv = 0.0
        prev_dt = self._adjustedFixedDates[start_index - 1]
        df_discount = 1.0
        if len(self._adjustedFixedDates) == 1:
            return 0.0

        for next_dt in self._adjustedFixedDates[start_index:]:
            alpha = day_counter.year_frac(prev_dt, next_dt)[0]
            df_discount = discount_curve.df(next_dt) / self._dfValuationDate
            flow = self._fixed_coupon * alpha * self._notional
            flowPV = flow * df_discount
            pv += flowPV
            prev_dt = next_dt

            self._fixed_year_fracs.append(alpha)
            self._fixedFlows.append(flow)
            self._fixedDfs.append(df_discount)
            self._fixedFlowPVs.append(flow * df_discount)
            self._fixed_total_pv.append(pv)

        flow = principal * self._notional
        pv = pv + flow * df_discount
        self._fixedFlowPVs[-1] += flow * df_discount
        self._fixedFlows[-1] += flow
        self._fixed_total_pv[-1] = pv
        return pv

##########################################################################

    def _calc_fixed_leg_flows(self):

        self._fixed_year_fracs = []
        self._fixedFlows = []

        day_counter = FinDayCount(self._fixed_day_count_type)

        """ Now PV fixed leg flows. """
        prev_dt = self._adjustedFixedDates[0]

        for next_dt in self._adjustedFixedDates[1:]:
            alpha = day_counter.year_frac(prev_dt, next_dt)[0]
            flow = self._fixed_coupon * alpha * self._notional
            prev_dt = next_dt
            self._fixed_year_fracs.append(alpha)
            self._fixedFlows.append(flow)

##########################################################################

    def cash_settled_pv01(self,
                          valuation_date,
                          flat_swap_rate,
                          frequency_type):
        """ Calculate the forward value of an annuity of a forward starting
        swap using a single flat discount rate equal to the swap rate. This is
        used in the pricing of a cash-settled swaption in the IborSwaption
        class. This method does not affect the standard valuation methods."""

        m = FinFrequency(frequency_type)

        if m == 0:
            raise FinError("Frequency cannot be zero.")

        """ The swap may have started in the past but we can only value
        payments that have occurred after the valuation date. """
        start_index = 0
        while self._adjustedFixedDates[start_index] < valuation_date:
            start_index += 1

        """ If the swap has yet to settle then we do not include the
        start date of the swap as a coupon payment date. """
        if valuation_date <= self._effective_date:
            start_index = 1

        """ Now PV fixed leg flows. """
        flatPV01 = 0.0
        df = 1.0
        alpha = 1.0 / m

        for _ in self._adjustedFixedDates[start_index:]:
            df = df / (1.0 + alpha * flat_swap_rate)
            flatPV01 += df * alpha

        return flatPV01

##########################################################################

    def float_leg_value(self,
                        valuation_date,  # This should be the settlement date
                        discount_curve,
                        index_curve,
                        first_fixing_rate=None,
                        principal=0.0):
        """ Value the floating leg with payments from an index curve and
        discounting based on a supplied discount curve. The valuation date can
        be the today date. In this case the price of the floating leg will not
        be par (assuming we added on a principal repayment). This is only the
        case if we set the valuation date to be the swap's actual settlement
        date. """

        self._valuation_date = valuation_date
        self._float_year_fracs = []
        self._floatFlows = []
        self._floatRates = []
        self._floatDfs = []
        self._floatFlowPVs = []
        self._floatTotalPV = []
        self._firstFixingRate = first_fixing_rate

        basis = FinDayCount(self._float_day_count_type)

        """ The swap may have started in the past but we can only value
        payments that have occurred after the start date. """
        start_index = 0
        while self._adjustedFloatDates[start_index] < valuation_date:
            start_index += 1

        """ If the swap has yet to settle then we do not include the
        start date of the swap as a coupon payment date. """
        if valuation_date <= self._effective_date:
            start_index = 1

        self._floatStartIndex = start_index

        # Forward price to settlement date (if valuation is settlement date)
        self._dfValuationDate = discount_curve.df(valuation_date)

        """ The first floating payment is usually already fixed so is
        not implied by the index curve. """
        prev_dt = self._adjustedFloatDates[start_index - 1]
        next_dt = self._adjustedFloatDates[start_index]
        alpha = basis.year_frac(prev_dt, next_dt)[0]
        # Cannot be pcd as has past
        df1_index = index_curve.df(self._effective_date)
        df2_index = index_curve.df(next_dt)

        float_rate = 0.0

        if self._firstFixingRate is None:
            fwd_rate = (df1_index / df2_index - 1.0) / alpha
            flow = (fwd_rate + self._float_spread) * alpha * self._notional
            float_rate = fwd_rate
        else:
            flow = self._firstFixingRate * alpha * self._notional
            float_rate = self._firstFixingRate

        # All discounting is done forward to the valuation date
        df_discount = discount_curve.df(next_dt) / self._dfValuationDate

        pv = flow * df_discount

        self._float_year_fracs.append(alpha)
        self._floatFlows.append(flow)
        self._floatRates.append(float_rate)
        self._floatDfs.append(df_discount)
        self._floatFlowPVs.append(flow * df_discount)
        self._floatTotalPV.append(pv)

        prev_dt = next_dt
        df1_index = index_curve.df(prev_dt)

        for next_dt in self._adjustedFloatDates[start_index + 1:]:
            alpha = basis.year_frac(prev_dt, next_dt)[0]
            df2_index = index_curve.df(next_dt)
            # The accrual factors cancel
            fwd_rate = (df1_index / df2_index - 1.0) / alpha
            flow = (fwd_rate + self._float_spread) * alpha * self._notional

            # All discounting is done forward to the valuation date
            df_discount = discount_curve.df(next_dt) / self._dfValuationDate

            pv += flow * df_discount
            df1_index = df2_index
            prev_dt = next_dt

            self._floatFlows.append(flow)
            self._float_year_fracs.append(alpha)
            self._floatRates.append(fwd_rate)
            self._floatDfs.append(df_discount)
            self._floatFlowPVs.append(flow * df_discount)
            self._floatTotalPV.append(pv)

        flow = principal * self._notional
        pv = pv + flow * df_discount
        self._floatFlows[-1] += flow
        self._floatFlowPVs[-1] += flow * df_discount
        self._floatTotalPV[-1] = pv

        return pv

##########################################################################

    def print_fixed_leg_pv(self):
        """ Prints the fixed leg dates, accrual factors, discount factors,
        cash amounts, their present value and their cumulative PV using the
        last valuation performed. """

        print("START DATE:", self._effective_date)
        print("MATURITY DATE:", self._maturity_date)
        print("COUPON (%):", self._fixed_coupon * 100)
        print("FIXED LEG FREQUENCY:", str(self._fixed_frequency_type))
        print("FIXED LEG DAY COUNT:", str(self._fixed_day_count_type))
        print("VALUATION DATE", self._valuation_date)

        if len(self._fixedFlows) == 0:
            print("Fixed Flows not calculated.")
            return

        header = "PAYMENT_DATE     YEAR_FRAC        FLOW         DF"
        header += "         DF*FLOW       CUM_PV"
        print(header)

        if self._fixedStartIndex is None:
            raise FinError("Need to value swap before calling this function.")

        start_index = self._fixedStartIndex

        # By definition the discount factor is 1.0 on the valuation date
        print("%15s %10s %12s %12.8f %12s %12s" %
              (self._valuation_date,
               "-",
               "-",
               1.0,
               "-",
               "-"))

        iFlow = 0
        for payment_date in self._adjustedFixedDates[start_index:]:
            print("%15s %10.7f %12.2f %12.8f %12.2f %12.2f" %
                  (payment_date,
                   self._fixed_year_fracs[iFlow],
                   self._fixedFlows[iFlow],
                   self._fixedDfs[iFlow],
                   self._fixedFlowPVs[iFlow],
                   self._fixed_total_pv[iFlow]))

            iFlow += 1

##########################################################################

    def print_fixed_leg_flows(self):
        """ Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows. """

        print("START DATE:", self._effective_date)
        print("MATURITY DATE:", self._maturity_date)
        print("COUPON (%):", self._fixed_coupon * 100)
        print("FIXED LEG FREQUENCY:", str(self._fixed_frequency_type))
        print("FIXED LEG DAY COUNT:", str(self._fixed_day_count_type))

        if len(self._fixedFlows) == 0:
            print("Fixed Flows not calculated.")
            return

        header = "PAYMENT_DATE     YEAR_FRAC        FLOW"
        print(header)

        start_index = 1

        iFlow = 0
        for payment_date in self._adjustedFixedDates[start_index:]:
            print("%15s %12.8f %12.2f" %
                  (payment_date,
                   self._fixed_year_fracs[iFlow],
                   self._fixedFlows[iFlow]))

            iFlow += 1

##########################################################################

    def print_float_leg_pv(self):
        """ Prints the floating leg dates, accrual factors, discount factors,
        forward libor rates, implied cash amounts, their present value and
        their cumulative PV using the last valuation performed. """

        print("START DATE:", self._effective_date)
        print("MATURITY DATE:", self._maturity_date)
        print("SPREAD COUPON (%):", self._float_spread * 100)
        print("FLOAT LEG FREQUENCY:", str(self._float_frequency_type))
        print("FLOAT LEG DAY COUNT:", str(self._float_day_count_type))
        print("VALUATION DATE", self._valuation_date)

        if len(self._floatFlows) == 0:
            print("Floating Flows not calculated.")
            return

        if self._firstFixingRate is None:
            print("         *** FIRST FLOATING RATE PAYMENT IS IMPLIED ***")

        header = "PAYMENT_DATE     YEAR_FRAC    RATE(%)       FLOW         DF"
        header += "         DF*FLOW       CUM_PV"
        print(header)

        start_index = self._floatStartIndex

        # By definition the discount factor is 1.0 on the valuation date

        print("%15s %10s %10s %12s %12.8f %12s %12s" %
              (self._valuation_date,
               "-",
               "-",
               "-",
               1.0,
               "-",
               "-"))

        iFlow = 0
        for payment_date in self._adjustedFloatDates[start_index:]:
            print("%15s %10.7f %10.5f %12.2f %12.8f %12.2f %12.2f" %
                  (payment_date,
                   self._float_year_fracs[iFlow],
                   self._floatRates[iFlow]*100.0,
                   self._floatFlows[iFlow],
                   self._floatDfs[iFlow],
                   self._floatFlowPVs[iFlow],
                   self._floatTotalPV[iFlow]))

            iFlow += 1

##########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("START DATE", self._effective_date)
        s += label_to_string("TERMINATION DATE", self._termination_date)
        s += label_to_string("MATURITY DATE", self._maturity_date)
        s += label_to_string("NOTIONAL", self._notional)
        s += label_to_string("SWAP TYPE", self._swapType)
        s += label_to_string("FIXED COUPON", self._fixed_coupon)
        s += label_to_string("FLOAT SPREAD", self._float_spread)
        s += label_to_string("FIXED FREQUENCY", self._fixed_frequency_type)
        s += label_to_string("FLOAT FREQUENCY", self._float_frequency_type)
        s += label_to_string("FIXED DAY COUNT", self._fixed_day_count_type)
        s += label_to_string("FLOAT DAY COUNT", self._float_day_count_type)
        s += label_to_string("CALENDAR", self._calendar_type)
        s += label_to_string("BUS DAY ADJUST", self._bus_day_adjust_type)
        s += label_to_string("DATE GEN TYPE", self._date_gen_rule_type)
        return s

###############################################################################

    def _print(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        print(self)

###############################################################################
