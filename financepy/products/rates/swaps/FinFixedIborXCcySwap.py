##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...utils.FinError import FinError
from ...utils.Date import Date
from ...utils.FinGlobalVariables import g_small
from ...utils.FinDayCount import FinDayCount, DayCountTypes
from ...utils.FinFrequency import FrequencyTypes, FinFrequency
from ...utils.FinCalendar import CalendarTypes,  DateGenRuleTypes
from ...utils.FinCalendar import Calendar, BusDayAdjustTypes
from ...utils.FinSchedule import FinSchedule
from ...utils.FinHelperFunctions import label_to_string, check_argument_types
from ...utils.FinMath import ONE_MILLION
from ...utils.FinGlobalTypes import SwapTypes

##########################################################################


class FinFixedIborXCcySwap:
    """ Class for managing a cross currency swap contract. This is a contract
    in which a fixed payment leg in one currency is exchanged for floating
    payments in a second currency. There is an exchange of par at maturity.
    The contract is entered into at zero initial cost and it lasts from a start
    date to a specified maturity date.

    The value of the contract is the NPV of the two cpn streams. Discounting
    is done on a supplied discount discount (one for each leg) which is separate
    from the curve from which the implied index rates are extracted. """

    def __init__(self,
                 effective_dt: Date,  # Date interest starts to accrue
                 term_dt_or_tenor: (Date, str),  # Date contract ends
                 fixed_leg_type: SwapTypes,
                 fixed_cpn: float,  # Fixed cpn (annualised)
                 fixed_freq_type: FrequencyTypes,
                 fixed_dc_type: DayCountTypes,
                 float_spread: float = 0.0,
                 float_freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 float_dc_type: DayCountTypes = DayCountTypes.THIRTY_E_360,
                 notional: float = ONE_MILLION,
                 cal_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):
        """ Create an interest rate swap contract giving the contract start
        date, its maturity, fixed cpn, fixed leg frequency, fixed leg day
        count convention and notional. The floating leg parameters have default
        values that can be overwritten if needed. The start date is contractual
        and is the same as the settlement date for a new swap. It is the date
        on which interest starts to accrue. The end of the contract is the
        termination date. This is not adjusted for business days. The adjusted
        termination date is called the maturity date. This is calculated. """

        check_argument_types(self.__init__, locals())

        if isinstance(term_dt_or_tenor, Date):
            self._termination_dt = term_dt_or_tenor
        else:
            self._termination_dt = effective_dt.add_tenor(
                term_dt_or_tenor)

        calendar = Calendar(cal_type)
        self.maturity_dt = calendar.adjust(self._termination_dt,
                                            bd_type)

        if effective_dt > self.maturity_dt:
            raise FinError("Start date after maturity date")

        self.effective_dt = effective_dt
        self.notional = notional

        self._fixed_cpn = fixed_cpn
        self._float_spread = float_spread

        self._fixed_freq_type = fixed_freq_type
        self._float_freq_type = float_freq_type

        self._fixed_dc_type = fixed_dc_type
        self._float_dc_type = float_dc_type

        self._fixed_leg_type = fixed_leg_type

        self._cal_type = cal_type
        self._bd_type = bd_type
        self._dg_type = dg_type

        # These are generated immediately as they are for the entire
        # life of the swap. Given a valuation date we can determine
        # which cash flows are in the future and value the swap
        self._generate_fixed_leg_payment_dts()
        self._generate_float_leg_payment_dts()

        self._adjustedMaturityDate = self._adjusted_fixed_dts[-1]

        # Need to know latest payment date for bootstrap - DO I NEED THIS ??!
        self._lastPaymentDate = self.maturity_dt
        if self._adjusted_fixed_dts[-1] > self._lastPaymentDate:
            self._lastPaymentDate = self._adjusted_fixed_dts[-1]

        if self._adjustedFloatDates[-1] > self._lastPaymentDate:
            self._lastPaymentDate = self._adjustedFloatDates[-1]

        # NOT TO BE PRINTED
        self._floatYearFracs = []
        self._float_flows = []
        self._float_flow_pvs = []
        self._float_dfs = []

        self._fixedYearFracs = []
        self._fixed_flows = []
        self._fixedDfs = []
        self._fixedFlowPVs = []

        self._first_fixing_rate = None
        self.value_dt = None
        self._fixedStartIndex = None

        self._calc_fixed_leg_flows()

##########################################################################

    def value(self,
              value_dt,
              discount_curve,
              index_curve,
              first_fixing_rate=None,
              principal=0.0):
        """ Value the interest rate swap on a value date given a single Libor
        discount curve. """

        fixed_leg_value = self.fixed_leg_value(value_dt,
                                               discount_curve,
                                               principal)

        float_leg_value = self.float_leg_value(value_dt,
                                               discount_curve,
                                               index_curve,
                                               first_fixing_rate,
                                               principal)

        value = fixed_leg_value - float_leg_value

        if self._fixed_leg_type == SwapTypes.PAY:
            value = value * (-1.0)

        return value

##########################################################################

    def _generate_fixed_leg_payment_dts(self):
        """ Generate the fixed leg payment dates all the way back to
        the start date of the swap which may precede the valuation date"""
        self._adjusted_fixed_dts = FinSchedule(
            self.effective_dt,
            self._termination_dt,
            self._fixed_freq_type,
            self._cal_type,
            self._bd_type,
            self._dg_type).generate()

##########################################################################

    def _generate_float_leg_payment_dts(self):
        """ Generate the floating leg payment dates all the way back to
        the start date of the swap which may precede the valuation date"""
        self._adjustedFloatDates = FinSchedule(
            self.effective_dt,
            self._termination_dt,
            self._float_freq_type,
            self._cal_type,
            self._bd_type,
            self._dg_type).generate()

##########################################################################

    def fixed_dts(self):
        """ return a vector of the fixed leg payment dates """
        if self._adjusted_fixed_dts is None:
            raise FinError("Fixed dates have not been generated")

        return self._adjusted_fixed_dts[1:]

##########################################################################

    def float_dts(self):
        """ return a vector of the fixed leg payment dates """
        if self._adjustedFloatDates is None:
            raise FinError("Float dates have not been generated")

        return self._adjustedFloatDates[1:]

##########################################################################

    def pv01(self, value_dt, discount_curve):
        """ Calculate the value of 1 basis point cpn on the fixed leg. """

        pv = self.fixed_leg_value(value_dt, discount_curve)
        pv01 = pv / self._fixed_cpn / self.notional
        return pv01

##########################################################################

    def swap_rate(self, value_dt, discount_curve):
        """ Calculate the fixed leg cpn that makes the swap worth zero.
        If the valuation date is before the swap payments start then this
        is the forward swap rate as it starts in the future. The swap rate
        is then a forward swap rate and so we use a forward discount
        factor. If the swap fixed leg has begun then we have a spot
        starting swap. """

        pv01 = self.pv01(value_dt, discount_curve)

        if value_dt < self.effective_dt:
            df0 = discount_curve.df(self.effective_dt)
        else:
            df0 = discount_curve.df(value_dt)

        dfT = discount_curve.df(self.maturity_dt)

        if abs(pv01) < g_small:
            raise FinError("PV01 is zero. Cannot compute swap rate.")

        cpn = (df0 - dfT) / pv01
        return cpn

##########################################################################

    def fixed_leg_value(self, value_dt, discount_curve, principal=0.0):

        self.value_dt = value_dt

        self._fixedYearFracs = []
        self._fixed_flows = []
        self._fixedDfs = []
        self._fixedFlowPVs = []
        self._fixed_total_pv = []

        day_counter = FinDayCount(self._fixed_dc_type)

        """ The swap may have started in the past but we can only value
        payments that have occurred after the valuation date. """
        start_index = 0
        while self._adjusted_fixed_dts[start_index] < value_dt:
            start_index += 1

        """ If the swap has yet to settle then we do not include the
        start date of the swap as a cpn payment date. """
        if value_dt <= self.effective_dt:
            start_index = 1

        self._fixedStartIndex = start_index

        """ Now PV fixed leg flows. """
        self._df_value_dt = discount_curve.df(value_dt)

        pv = 0.0
        prev_dt = self._adjusted_fixed_dts[start_index - 1]
        df_discount = 1.0
        if len(self._adjusted_fixed_dts) == 1:
            return 0.0

        for next_dt in self._adjusted_fixed_dts[start_index:]:
            alpha = day_counter.year_frac(prev_dt, next_dt)[0]
            df_discount = discount_curve.df(next_dt) / self._df_value_dt
            flow = self._fixed_cpn * alpha * self.notional
            flowPV = flow * df_discount
            pv += flowPV
            prev_dt = next_dt

            self._fixedYearFracs.append(alpha)
            self._fixed_flows.append(flow)
            self._fixedDfs.append(df_discount)
            self._fixedFlowPVs.append(flow * df_discount)
            self._fixed_total_pv.append(pv)

        flow = principal * self.notional
        pv = pv + flow * df_discount
        self._fixedFlowPVs[-1] += flow * df_discount
        self._fixed_flows[-1] += flow
        self._fixed_total_pv[-1] = pv
        return pv

##########################################################################

    def _calc_fixed_leg_flows(self):

        self._fixedYearFracs = []
        self._fixed_flows = []

        day_counter = FinDayCount(self._fixed_dc_type)

        """ Now PV fixed leg flows. """
        prev_dt = self._adjusted_fixed_dts[0]

        for next_dt in self._adjusted_fixed_dts[1:]:
            alpha = day_counter.year_frac(prev_dt, next_dt)[0]
            flow = self._fixed_cpn * alpha * self.notional
            prev_dt = next_dt
            self._fixedYearFracs.append(alpha)
            self._fixed_flows.append(flow)

##########################################################################

    def cash_settled_pv01(self,
                          value_dt,
                          flatSwapRate,
                          frequencyType):
        """ Calculate the forward value of an annuity of a forward starting
        swap using a single flat discount rate equal to the swap rate. This is
        used in the pricing of a cash-settled swaption in the IborSwaption
        class. This method does not affect the standard valuation methods."""

        m = FinFrequency(frequencyType)

        if m == 0:
            raise FinError("Frequency cannot be zero.")

        """ The swap may have started in the past but we can only value
        payments that have occurred after the valuation date. """
        start_index = 0
        while self._adjusted_fixed_dts[start_index] < value_dt:
            start_index += 1

        """ If the swap has yet to settle then we do not include the
        start date of the swap as a cpn payment date. """
        if value_dt <= self.effective_dt:
            start_index = 1

        """ Now PV fixed leg flows. """
        flat_pv01 = 0.0
        df = 1.0
        alpha = 1.0 / m

        for _ in self._adjusted_fixed_dts[start_index:]:
            df = df / (1.0 + alpha * flatSwapRate)
            flat_pv01 += df * alpha

        return flat_pv01

##########################################################################

    def float_leg_value(self,
                        value_dt,  # This should be the settlement date
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

        self.value_dt = value_dt
        self._floatYearFracs = []
        self._float_flows = []
        self._floatRates = []
        self._float_dfs = []
        self._float_flow_pvs = []
        self._floatTotalPV = []
        self._first_fixing_rate = first_fixing_rate

        basis = FinDayCount(self._float_dc_type)

        """ The swap may have started in the past but we can only value
        payments that have occurred after the start date. """
        start_index = 0
        while self._adjustedFloatDates[start_index] < value_dt:
            start_index += 1

        """ If the swap has yet to settle then we do not include the
        start date of the swap as a cpn payment date. """
        if value_dt <= self.effective_dt:
            start_index = 1

        self._floatStartIndex = start_index

        # Forward price to settlement date (if valuation is settlement date)
        self._df_value_dt = discount_curve.df(value_dt)

        """ The first floating payment is usually already fixed so is
        not implied by the index curve. """
        prev_dt = self._adjustedFloatDates[start_index - 1]
        next_dt = self._adjustedFloatDates[start_index]
        alpha = basis.year_frac(prev_dt, next_dt)[0]
        # Cannot be pcd as has past
        df1_index = index_curve.df(self.effective_dt)
        df2_index = index_curve.df(next_dt)

        floatRate = 0.0

        if self._first_fixing_rate is None:
            fwd_rate = (df1_index / df2_index - 1.0) / alpha
            flow = (fwd_rate + self._float_spread) * alpha * self.notional
            floatRate = fwd_rate
        else:
            flow = self._first_fixing_rate * alpha * self.notional
            floatRate = self._first_fixing_rate

        # All discounting is done forward to the valuation date
        df_discount = discount_curve.df(next_dt) / self._df_value_dt

        pv = flow * df_discount

        self._floatYearFracs.append(alpha)
        self._float_flows.append(flow)
        self._floatRates.append(floatRate)
        self._float_dfs.append(df_discount)
        self._float_flow_pvs.append(flow * df_discount)
        self._floatTotalPV.append(pv)

        prev_dt = next_dt
        df1_index = index_curve.df(prev_dt)

        for next_dt in self._adjustedFloatDates[start_index + 1:]:
            alpha = basis.year_frac(prev_dt, next_dt)[0]
            df2_index = index_curve.df(next_dt)
            # The accrual factors cancel
            fwd_rate = (df1_index / df2_index - 1.0) / alpha
            flow = (fwd_rate + self._float_spread) * alpha * self.notional

            # All discounting is done forward to the valuation date
            df_discount = discount_curve.df(next_dt) / self._df_value_dt

            pv += flow * df_discount
            df1_index = df2_index
            prev_dt = next_dt

            self._float_flows.append(flow)
            self._floatYearFracs.append(alpha)
            self._floatRates.append(fwd_rate)
            self._float_dfs.append(df_discount)
            self._float_flow_pvs.append(flow * df_discount)
            self._floatTotalPV.append(pv)

        flow = principal * self.notional
        pv = pv + flow * df_discount
        self._float_flows[-1] += flow
        self._float_flow_pvs[-1] += flow * df_discount
        self._floatTotalPV[-1] = pv

        return pv

##########################################################################

    def print_fixed_leg_pv(self):
        """ Prints the fixed leg dates, accrual factors, discount factors,
        cash amounts, their present value and their cumulative PV using the
        last valuation performed. """

        print("START DATE:", self.effective_dt)
        print("MATURITY DATE:", self.maturity_dt)
        print("cpn (%):", self._fixed_cpn * 100)
        print("FIXED LEG FREQUENCY:", str(self._fixed_freq_type))
        print("FIXED LEG DAY COUNT:", str(self._fixed_dc_type))
        print("VALUATION DATE", self.value_dt)

        if len(self._fixed_flows) == 0:
            print("Fixed Flows not calculated.")
            return

        header = "PAYMENT_dt     YEAR_FRAC        FLOW         DF"
        header += "         DF*FLOW       CUM_PV"
        print(header)

        if self._fixedStartIndex is None:
            raise FinError("Need to value swap before calling this function.")

        start_index = self._fixedStartIndex

        # By definition the discount factor is 1.0 on the valuation date
        print("%15s %10s %12s %12.8f %12s %12s" %
              (self.value_dt,
               "-",
               "-",
               1.0,
               "-",
               "-"))

        i_flow = 0
        for payment_dt in self._adjusted_fixed_dts[start_index:]:
            print("%15s %10.7f %12.2f %12.8f %12.2f %12.2f" %
                  (payment_dt,
                   self._fixedYearFracs[i_flow],
                   self._fixed_flows[i_flow],
                   self._fixedDfs[i_flow],
                   self._fixedFlowPVs[i_flow],
                   self._fixed_total_pv[i_flow]))

            i_flow += 1

##########################################################################

    def print_fixed_leg_flows(self):
        """ Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows. """

        print("START DATE:", self.effective_dt)
        print("MATURITY DATE:", self.maturity_dt)
        print("cpn (%):", self._fixed_cpn * 100)
        print("FIXED LEG FREQUENCY:", str(self._fixed_freq_type))
        print("FIXED LEG DAY COUNT:", str(self._fixed_dc_type))

        if len(self._fixed_flows) == 0:
            print("Fixed Flows not calculated.")
            return

        header = "PAYMENT_dt     YEAR_FRAC        FLOW"
        print(header)

        start_index = 1

        i_flow = 0
        for payment_dt in self._adjusted_fixed_dts[start_index:]:
            print("%15s %12.8f %12.2f" %
                  (payment_dt,
                   self._fixedYearFracs[i_flow],
                   self._fixed_flows[i_flow]))

            i_flow += 1

##########################################################################

    def print_float_leg_pv(self):
        """ Prints the floating leg dates, accrual factors, discount factors,
        forward libor rates, implied cash amounts, their present value and
        their cumulative PV using the last valuation performed. """

        print("START DATE:", self.effective_dt)
        print("MATURITY DATE:", self.maturity_dt)
        print("SPREAD cpn (%):", self._float_spread * 100)
        print("FLOAT LEG FREQUENCY:", str(self._float_freq_type))
        print("FLOAT LEG DAY COUNT:", str(self._float_dc_type))
        print("VALUATION DATE", self.value_dt)

        if len(self._float_flows) == 0:
            print("Floating Flows not calculated.")
            return

        if self._first_fixing_rate is None:
            print("         *** FIRST FLOATING RATE PAYMENT IS IMPLIED ***")

        header = "PAYMENT_dt     YEAR_FRAC    RATE(%)       FLOW         DF"
        header += "         DF*FLOW       CUM_PV"
        print(header)

        start_index = self._floatStartIndex

        # By definition the discount factor is 1.0 on the valuation date

        print("%15s %10s %10s %12s %12.8f %12s %12s" %
              (self.value_dt,
               "-",
               "-",
               "-",
               1.0,
               "-",
               "-"))

        i_flow = 0
        for payment_dt in self._adjustedFloatDates[start_index:]:
            print("%15s %10.7f %10.5f %12.2f %12.8f %12.2f %12.2f" %
                  (payment_dt,
                   self._floatYearFracs[i_flow],
                   self._floatRates[i_flow]*100.0,
                   self._float_flows[i_flow],
                   self._float_dfs[i_flow],
                   self._float_flow_pvs[i_flow],
                   self._floatTotalPV[i_flow]))

            i_flow += 1

##########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("START DATE", self.effective_dt)
        s += label_to_string("TERMINATION DATE", self._termination_dt)
        s += label_to_string("MATURITY DATE", self.maturity_dt)
        s += label_to_string("NOTIONAL", self.notional)
        s += label_to_string("SWAP TYPE", self._swap_type)
        s += label_to_string("FIXED cpn", self._fixed_cpn)
        s += label_to_string("FLOAT SPREAD", self._float_spread)
        s += label_to_string("FIXED FREQUENCY", self._fixed_freq_type)
        s += label_to_string("FLOAT FREQUENCY", self._float_freq_type)
        s += label_to_string("FIXED DAY COUNT", self._fixed_dc_type)
        s += label_to_string("FLOAT DAY COUNT", self._float_dc_type)
        s += label_to_string("CALENDAR", self._cal_type)
        s += label_to_string("BUS DAY ADJUST", self._bd_type)
        s += label_to_string("DATE GEN TYPE", self._dg_type)
        return s

###############################################################################

    def _print(self):
        """ Print a list of the unadjusted cpn payment dates used in
        analytic calculations for the bond. """
        print(self)

###############################################################################
