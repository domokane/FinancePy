##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...utils.FinError import FinError
from ...utils.Date import Date
from ...utils.FinDayCount import FinDayCount, DayCountTypes
from ...utils.FinFrequency import FrequencyTypes, FinFrequency
from ...utils.FinCalendar import CalendarTypes, DateGenRuleTypes
from ...utils.FinCalendar import Calendar, BusDayAdjustTypes
from ...utils.FinSchedule import FinSchedule
from ...utils.FinHelperFunctions import label_to_string, check_argument_types
from ...utils.FinMath import ONE_MILLION

##########################################################################


class IborIborSwap:
    """ Class for managing an interest rate basis swap contract. This is a
    contract in which a floating leg with one LIBOR tenor is exchanged for a
    floating leg payment in a different LIBOR tenor. There is no exchange of
    par. The contract is entered into at zero initial cost. The contract lasts
    from a start date to a specified maturity date.

    The value of the contract is the NPV of the two cpn streams. Discounting
    is done on a supplied discount curve which is separate from the discount from
    which the implied index rates are extracted. """

    def __init__(self,
                 effective_dt: Date,  # Date interest starts to accrue
                 term_dt_or_tenor: (Date, str),  # Date contract ends
                 payFreqType: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 payDayCountType: DayCountTypes = DayCountTypes.THIRTY_E_360,
                 recFreqType: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 recDayCountType: DayCountTypes = DayCountTypes.THIRTY_E_360,
                 basisSwapSpread: float = 0.0,
                 notional: float = ONE_MILLION,
                 cal_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):
        """ Create an Ibor basis swap contract giving the contract start
        date, its maturity, frequency and day counts on the two floating
        legs and notional. The floating leg parameters have default
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

        calendar = CalendarTypes(cal_type)
        self.maturity_dt = calendar.adjust(self._termination_dt,
                                            bd_type)

        if effective_dt > self.maturity_dt:
            raise FinError("Start date after maturity date")

        self.effective_dt = effective_dt
        self.notional = notional
        self._basisSwapSpread = basisSwapSpread
        self._payFreqType = payFreqType
        self._recFreqType = recFreqType
        self._payDayCountType = payDayCountType
        self._recDayCountType = recDayCountType
        self._cal_type = cal_type
        self._bd_type = bd_type
        self._dg_type = dg_type

        self._payFloatDates = self._generateFloatLegDates(self._payFreqType)
        self._recFloatDates = self._generateFloatLegDates(self._recFreqType)

        self._adjustedMaturityDate = self._adjusted_fixed_dts[-1]

        self._payFloatYearFracs = []
        self._recFloatYearFracs = []

        self._payfloat_flows = []
        self._recfloat_flows = []

        self._pay_float_flow_pvs = []
        self._recfloat_flow_pvs = []

        self._payfirst_fixing_rate = None
        self._recfirst_fixing_rate = None

        self.value_dt = None

##########################################################################

    def _generate_pay_float_leg_payment_dts(self, freq_type):
        """ Generate the floating leg payment dates all the way back to
        the start date of the swap which may precede the valuation date"""

        floatDates = FinSchedule(self.effective_dt,
                                 self._termination_dt,
                                 freq_type,
                                 self._cal_type,
                                 self._bd_type,
                                 self._dg_type).generate()

        return floatDates

##########################################################################

    def value(self,
              value_dt,
              discount_curve,
              payIndexCurve,
              recIndexCurve,
              payfirst_fixing_rate=None,
              recfirst_fixing_rate=None,
              principal=0.0):
        """ Value the LIBOR basis swap on a value date given a single Ibor
        discount curve and each of the index discount for the two floating legs
        of the swap. """

        payFloatLegValue = self.float_leg_value(value_dt,
                                                discount_curve,
                                                payIndexCurve,
                                                payfirst_fixing_rate,
                                                principal)

        recFloatLegValue = self.float_leg_value(value_dt,
                                                discount_curve,
                                                recIndexCurve,
                                                recfirst_fixing_rate,
                                                principal)

        value = recFloatLegValue - payFloatLegValue
        return value

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
        self._floatDfs = []
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
        self._floatDfs.append(df_discount)
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
            self._floatDfs.append(df_discount)
            self._float_flow_pvs.append(flow * df_discount)
            self._floatTotalPV.append(pv)

        flow = principal * self.notional
        pv = pv + flow * df_discount
        self._float_flows[-1] += flow
        self._float_flow_pvs[-1] += flow * df_discount
        self._floatTotalPV[-1] = pv

        return pv

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
                   self._floatDfs[i_flow],
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
        s += label_to_string("FIXED COUPON", self._fixed_cpn)
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
