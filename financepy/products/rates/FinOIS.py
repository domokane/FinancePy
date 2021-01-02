##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinCalendar import FinCalendarTypes,  FinDateGenRuleTypes
from ...finutils.FinCalendar import FinCalendar, FinBusDayAdjustTypes
from ...finutils.FinHelperFunctions import checkArgumentTypes
from ...finutils.FinMath import ONE_MILLION
from ...finutils.FinGlobalTypes import FinSwapTypes
from ...market.curves.FinDiscountCurve import FinDiscountCurve

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
    ''' Class for managing overnight index rate swaps (OIS) and Fed Fund swaps. 
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
    market OIS rates. '''

    def __init__(self,
                 effectiveDate: FinDate,  # Date interest starts to accrue
                 terminationDateOrTenor: (FinDate, str),  # Date contract ends
                 fixedLegType: FinSwapTypes,
                 fixedCoupon: float,  # Fixed coupon (annualised)
                 fixedFreqType: FinFrequencyTypes,
                 fixedDayCountType: FinDayCountTypes,
                 notional: float = ONE_MILLION,
                 paymentLag: int = 0, # Number of days after period payment occurs
                 floatSpread: float = 0.0,
                 floatFreqType: FinFrequencyTypes = FinFrequencyTypes.ANNUAL,
                 floatDayCountType: FinDayCountTypes = FinDayCountTypes.THIRTY_E_360,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD):
        ''' Create an overnight index swap contract giving the contract start
        date, its maturity, fixed coupon, fixed leg frequency, fixed leg day
        count convention and notional. The floating leg parameters have default
        values that can be overwritten if needed. The start date is contractual
        and is the same as the settlement date for a new swap. It is the date
        on which interest starts to accrue. The end of the contract is the
        termination date. This is not adjusted for business days. The adjusted
        termination date is called the maturity date. This is calculated. '''

        checkArgumentTypes(self.__init__, locals())

        if type(terminationDateOrTenor) == FinDate:
            self._terminationDate = terminationDateOrTenor
        else:
            self._terminationDate = effectiveDate.addTenor(terminationDateOrTenor)

        calendar = FinCalendar(calendarType)
        self._maturityDate = calendar.adjust(self._terminationDate,
                                             busDayAdjustType)

        if effectiveDate > self._maturityDate:
            raise FinError("Start date after maturity date")

        self._effectiveDate = effectiveDate

        floatLegType = FinSwapTypes.PAY
        if fixedLegType == FinSwapTypes.PAY:
            floatLegType = FinSwapTypes.RECEIVE

        principal = 0.0

        self._fixedLeg = FinFixedLeg(effectiveDate,
                                     self._terminationDate,
                                     fixedLegType,
                                     fixedCoupon,
                                     fixedFreqType,
                                     fixedDayCountType,
                                     notional,
                                     principal,
                                     paymentLag,
                                     calendarType,
                                     busDayAdjustType,
                                     dateGenRuleType)

        self._floatLeg = FinFloatLeg(effectiveDate,
                                     self._terminationDate,
                                     floatLegType,
                                     floatSpread,
                                     floatFreqType,
                                     floatDayCountType,
                                     notional,
                                     principal,
                                     paymentLag,
                                     calendarType,
                                     busDayAdjustType,
                                     dateGenRuleType)

###############################################################################

    def value(self,
              valuationDate: FinDate,
              oisCurve: FinDiscountCurve,
              firstFixingRate=None):
        ''' Value the interest rate swap on a value date given a single Ibor
        discount curve. '''

        fixedLegValue = self._fixedLeg.value(valuationDate,
                                             oisCurve)

        floatLegValue = self._floatLeg.value(valuationDate,
                                             oisCurve,
                                             oisCurve,
                                             firstFixingRate)

        value = fixedLegValue + floatLegValue
        return value
            
##########################################################################

    def pv01(self, valuationDate, discountCurve):
        ''' Calculate the value of 1 basis point coupon on the fixed leg. '''

        pv = self._fixedLeg.value(valuationDate, discountCurve)
        
        # Needs to be positive even if it is a payer leg
        pv = np.abs(pv)
        pv01 = pv / self._fixedLeg._coupon / self._fixedLeg._notional
        return pv01

##########################################################################

    def swapRate(self, valuationDate, oisCurve):
        ''' Calculate the fixed leg coupon that makes the swap worth zero.
        If the valuation date is before the swap payments start then this
        is the forward swap rate as it starts in the future. The swap rate
        is then a forward swap rate and so we use a forward discount
        factor. If the swap fixed leg has begun then we have a spot
        starting swap. '''

        pv01 = self.pv01(valuationDate, oisCurve)

        if valuationDate < self._effectiveDate:
            df0 = oisCurve.df(self._effectiveDate)
        else:
            df0 = oisCurve.df(valuationDate)

        floatLegPV = 0.0

        dfT = oisCurve.df(self._maturityDate)
        floatLegPV = (df0 - dfT) 
        floatLegPV /= self._fixedLeg._notional
        cpn = floatLegPV / pv01           
        return cpn

###############################################################################

    def printFixedLegPV(self):
        ''' Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows. '''

        self._fixedLeg.printValuation()

###############################################################################

    def printFloatLegPV(self):
        ''' Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows. '''

        self._floatLeg.printValuation()

###############################################################################

    def printFlows(self):
        ''' Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows. '''

        self._fixedLeg.printPayments()
        
        self._floatLeg.printPayments()

##########################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += self._fixedLeg.__repr__()
        s += "\n"
        s += self._floatLeg.__repr__()
        return s

###############################################################################

    def _print(self):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''
        print(self)

###############################################################################

