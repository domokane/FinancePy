##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...utils.FinError import FinError
from ...utils.date import Date
from ...utils.day_count import DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.calendar import CalendarTypes,  DateGenRuleTypes
from ...utils.calendar import Calendar, BusDayAdjustTypes
from ...utils.helper_functions import check_argument_types
from ...utils.fin_math import ONE_MILLION
from ...utils.FinGlobalTypes import FinSwapTypes
from ...market.curves.discount_curve import DiscountCurve

from .FinFixedLeg import FinFixedLeg
from .FinFloatLeg import FinFloatLeg

###############################################################################

from enum import Enum

class FinCompoundingTypes(Enum):
     COMPOUNDED = 1
     OVERNIGHT_COMPOUNDED_ANNUAL_RATE = 2
     AVERAGED = 3
     AVERAGED_DAILY = 4


###############################################################################

class FinOIS(object):
    """ Class for managing overnight index rate swaps (OIS) and Fed Fund swaps. 
    This is a contract in which a fixed payment leg is exchanged for a payment
    which pays the rolled-up overnight index rate (OIR). There is no exchange
    of par. The contract is entered into at zero initial cost.

    NOTE: This class is almost identical to FinIborSwap but will possibly
    deviate as distinctions between the two become clear to me. If not they 
    will be converged (or inherited) to avoid duplication.
    
    The contract lasts from a start date to a specified maturity date.
    The fixed coupon is the OIS fixed rate for the corresponding tenor which is
    set at contract initiation.

    The floating rate is not known fully until the end of each payment period.
    It's calculated at the contract maturity and is based on daily observations
    of the overnight index rate which are compounded according to a specific
    convention. Hence the OIS floating rate is determined by the history of the
    OIS rates.

    In its simplest form, there is just one fixed rate payment and one floating
    rate payment at contract maturity. However when the contract becomes longer
    than one year the floating and fixed payments become periodic, usually with
    annual exchanges of cash.

    The value of the contract is the NPV of the two coupon streams. Discounting
    is done on the OIS curve which is itself implied by the term structure of
    market OIS rates. """

    def __init__(self,
                 effective_date: Date,  # Date interest starts to accrue
                 termination_date_or_tenor: (Date, str),  # Date contract ends
                 fixed_legType: FinSwapTypes,
                 fixedCoupon: float,  # Fixed coupon (annualised)
                 fixedFreqType: FrequencyTypes,
                 fixedDayCountType: DayCountTypes,
                 notional: float = ONE_MILLION,
                 payment_lag: int = 0,  # Number of days after period payment occurs
                 floatSpread: float = 0.0,
                 floatFreqType: FrequencyTypes = FrequencyTypes.ANNUAL,
                 floatDayCountType: DayCountTypes = DayCountTypes.THIRTY_E_360,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):
        """ Create an overnight index swap contract giving the contract start
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
            self._termination_date = effective_date.addTenor(termination_date_or_tenor)

        calendar = Calendar(calendar_type)
        self._maturity_date = calendar.adjust(self._termination_date,
                                             bus_day_adjust_type)

        if effective_date > self._maturity_date:
            raise FinError("Start date after maturity date")

        self._effective_date = effective_date

        floatLegType = FinSwapTypes.PAY
        if fixed_legType == FinSwapTypes.PAY:
            floatLegType = FinSwapTypes.RECEIVE

        principal = 0.0

        self._fixed_leg = FinFixedLeg(effective_date,
                                     self._termination_date,
                                     fixed_legType,
                                     fixedCoupon,
                                     fixedFreqType,
                                     fixedDayCountType,
                                     notional,
                                     principal,
                                     payment_lag,
                                     calendar_type,
                                     bus_day_adjust_type,
                                     date_gen_rule_type)

        self._floatLeg = FinFloatLeg(effective_date,
                                     self._termination_date,
                                     floatLegType,
                                     floatSpread,
                                     floatFreqType,
                                     floatDayCountType,
                                     notional,
                                     principal,
                                     payment_lag,
                                     calendar_type,
                                     bus_day_adjust_type,
                                     date_gen_rule_type)

###############################################################################

    def value(self,
              valuation_date: Date,
              oisCurve: DiscountCurve,
              firstFixingRate=None):
        """ Value the interest rate swap on a value date given a single Ibor
        discount curve. """

        fixed_legValue = self._fixed_leg.value(valuation_date,
                                             oisCurve)

        floatLegValue = self._floatLeg.value(valuation_date,
                                             oisCurve,
                                             oisCurve,
                                             firstFixingRate)

        value = fixed_legValue + floatLegValue
        return value
            
##########################################################################

    def pv01(self, valuation_date, discount_curve):
        """ Calculate the value of 1 basis point coupon on the fixed leg. """

        pv = self._fixed_leg.value(valuation_date, discount_curve)
        
        # Needs to be positive even if it is a payer leg
        pv = np.abs(pv)
        pv01 = pv / self._fixed_leg._coupon / self._fixed_leg._notional
        return pv01

##########################################################################

    def swap_rate(self, valuation_date, oisCurve):
        """ Calculate the fixed leg coupon that makes the swap worth zero.
        If the valuation date is before the swap payments start then this
        is the forward swap rate as it starts in the future. The swap rate
        is then a forward swap rate and so we use a forward discount
        factor. If the swap fixed leg has begun then we have a spot
        starting swap. """

        pv01 = self.pv01(valuation_date, oisCurve)

        if valuation_date < self._effective_date:
            df0 = oisCurve.df(self._effective_date)
        else:
            df0 = oisCurve.df(valuation_date)

        floatLegPV = 0.0

        dfT = oisCurve.df(self._maturity_date)
        floatLegPV = (df0 - dfT) 
        floatLegPV /= self._fixed_leg._notional
        cpn = floatLegPV / pv01           
        return cpn

###############################################################################

    def printFixedLegPV(self):
        """ Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows. """

        self._fixed_leg.printValuation()

###############################################################################

    def printFloatLegPV(self):
        """ Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows. """

        self._floatLeg.printValuation()

###############################################################################

    def printFlows(self):
        """ Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows. """

        self._fixed_leg.printPayments()
        
        self._floatLeg.printPayments()

##########################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += self._fixed_leg.__repr__()
        s += "\n"
        s += self._floatLeg.__repr__()
        return s

###############################################################################

    def _print(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        print(self)

###############################################################################

