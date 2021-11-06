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


class IborIborSwap:
    """ Class for managing an interest rate basis swap contract. This is a
    contract in which a floating leg with one LIBOR tenor is exchanged for a 
    floating leg payment in a different LIBOR tenor. There is no exchange of
    par. The contract is entered into at zero initial cost. The contract lasts
    from a start date to a specified maturity date.

    The value of the contract is the NPV of the two coupon streams. Discounting
    is done on a supplied discount curve which is separate from the discount from
    which the implied index rates are extracted. """

    def __init__(self,
                 effective_date: Date,  # Date interest starts to accrue
                 termination_date_or_tenor: (Date, str),  # Date contract ends
                 payFreqType: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 payDayCountType: DayCountTypes = DayCountTypes.THIRTY_E_360,
                 recFreqType: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 recDayCountType: DayCountTypes = DayCountTypes.THIRTY_E_360,
                 basisSwapSpread: float = 0.0,
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

        self._effective_date = effective_date
        self._notional = notional
        self._basisSwapSpread = basisSwapSpread
        self._payFreqType = payFreqType
        self._recFreqType = recFreqType
        self._payDayCountType = payDayCountType
        self._recDayCountType = recDayCountType
        self._calendar_type = calendar_type
        self._bus_day_adjust_type = bus_day_adjust_type
        self._date_gen_rule_type = date_gen_rule_type

        self._payFloatDates = self._generateFloatLegDates(self._payFreqType)
        self._recFloatDates = self._generateFloatLegDates(self._recFreqType)

        self._adjustedMaturityDate = self._adjustedFixedDates[-1]

        self._payFloatYearFracs = []
        self._recFloatYearFracs = []

        self._payFloatFlows = []
        self._recFloatFlows = []

        self._payFloatFlowPVs = []
        self._recFloatFlowPVs = []

        self._payFirstFixingRate = None
        self._recFirstFixingRate = None

        self._valuation_date = None

##########################################################################

    def _generate_pay_float_leg_payment_dates(self, freq_type):
        """ Generate the floating leg payment dates all the way back to
        the start date of the swap which may precede the valuation date"""

        floatDates = FinSchedule(self._effective_date,
                                 self._termination_date,
                                 freq_type,
                                 self._calendar_type,
                                 self._bus_day_adjust_type,
                                 self._date_gen_rule_type)._generate()

        return floatDates

##########################################################################

    def value(self,
              valuation_date,
              discount_curve,
              payIndexCurve,
              recIndexCurve,
              payFirstFixingRate=None,
              recFirstFixingRate=None,
              principal=0.0):
        """ Value the LIBOR basis swap on a value date given a single Ibor
        discount curve and each of the index discount for the two floating legs
        of the swap. """

        payFloatLegValue = self.float_leg_value(valuation_date,
                                                discount_curve,
                                                payIndexCurve,
                                                payFirstFixingRate,
                                                principal)

        recFloatLegValue = self.float_leg_value(valuation_date,
                                                discount_curve,
                                                recIndexCurve,
                                                recFirstFixingRate,
                                                principal)

        value = recFloatLegValue - payFloatLegValue
        return value

##########################################################################

    def float_leg_value(self,
                        valuation_date,  # This should be the settlement date
                        discount_curve,
                        index_curve,
                        firstFixingRate=None,
                        principal=0.0):
        """ Value the floating leg with payments from an index curve and
        discounting based on a supplied discount curve. The valuation date can
        be the today date. In this case the price of the floating leg will not
        be par (assuming we added on a principal repayment). This is only the
        case if we set the valuation date to be the swap's actual settlement
        date. """

        self._valuation_date = valuation_date
        self._floatYearFracs = []
        self._floatFlows = []
        self._floatRates = []
        self._floatDfs = []
        self._floatFlowPVs = []
        self._floatTotalPV = []
        self._firstFixingRate = firstFixingRate

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

        floatRate = 0.0

        if self._firstFixingRate is None:
            fwd_rate = (df1_index / df2_index - 1.0) / alpha
            flow = (fwd_rate + self._float_spread) * alpha * self._notional
            floatRate = fwd_rate
        else:
            flow = self._firstFixingRate * alpha * self._notional
            floatRate = self._firstFixingRate

        # All discounting is done forward to the valuation date
        df_discount = discount_curve.df(next_dt) / self._dfValuationDate

        pv = flow * df_discount

        self._floatYearFracs.append(alpha)
        self._floatFlows.append(flow)
        self._floatRates.append(floatRate)
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
            self._floatYearFracs.append(alpha)
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
                   self._floatYearFracs[iFlow],
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
