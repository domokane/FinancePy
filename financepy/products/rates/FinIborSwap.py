##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinGlobalVariables import gSmall
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes, FinFrequency
from ...finutils.FinCalendar import FinCalendarTypes,  FinDateGenRuleTypes
from ...finutils.FinCalendar import FinCalendar, FinBusDayAdjustTypes
from ...finutils.FinHelperFunctions import checkArgumentTypes, labelToString
from ...finutils.FinMath import ONE_MILLION
from ...finutils.FinGlobalTypes import FinSwapTypes
from ...market.curves.FinDiscountCurve import FinDiscountCurve

from .FinFixedLeg import FinFixedLeg
from .FinFloatLeg import FinFloatLeg

##########################################################################


class FinIborSwap(object):
    ''' Class for managing a standard Fixed vs IBOR swap. This is a contract
    in which a fixed payment leg is exchanged for a series of floating rates
    payments linked to some IBOR index rate. There is no exchange of principal.
    The contract is entered into at zero initial cost. The contract lasts from
    a start date to a specified maturity date.

    The floating rate is not known fully until the end of the preceding payment
    period. It is set in advance and paid in arrears. 
    
    The value of the contract is the NPV of the two coupon streams. Discounting
    is done on a supplied discount curve which is separate from the curve from
    which the implied index rates are extracted. '''
    
    def __init__(self,
                 effectiveDate: FinDate,  # Date interest starts to accrue
                 terminationDateOrTenor: (FinDate, str),  # Date contract ends
                 fixedLegType: FinSwapTypes,
                 fixedCoupon: float,  # Fixed coupon (annualised)
                 fixedFreqType: FinFrequencyTypes,
                 fixedDayCountType: FinDayCountTypes,
                 notional: float = ONE_MILLION,
                 floatSpread: float = 0.0,
                 floatFreqType: FinFrequencyTypes = FinFrequencyTypes.QUARTERLY,
                 floatDayCountType: FinDayCountTypes = FinDayCountTypes.THIRTY_E_360,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD):
        ''' Create an interest rate swap contract giving the contract start
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
        
        paymentLag = 0
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
              discountCurve: FinDiscountCurve,
              indexCurve: FinDiscountCurve=None,
              firstFixingRate=None):
        ''' Value the interest rate swap on a value date given a single Ibor
        discount curve. '''

        if indexCurve is None:
            indexCurve = discountCurve

        fixedLegValue = self._fixedLeg.value(valuationDate,
                                             discountCurve)

        floatLegValue = self._floatLeg.value(valuationDate,
                                             discountCurve,
                                             indexCurve,
                                             firstFixingRate)

        value = fixedLegValue + floatLegValue
        return value

###############################################################################

    def pv01(self, valuationDate, discountCurve):
        ''' Calculate the value of 1 basis point coupon on the fixed leg. '''

        pv = self._fixedLeg.value(valuationDate, discountCurve)
        
        # Needs to be positive even if it is a payer leg
        pv = np.abs(pv)
        pv01 = pv / self._fixedLeg._coupon / self._fixedLeg._notional
        return pv01

###############################################################################

    def swapRate(self, 
                 valuationDate:FinDate,
                 discountCurve: FinDiscountCurve,
                 indexCurve: FinDiscountCurve = None,
                 firstFixing: float = None):
        ''' Calculate the fixed leg coupon that makes the swap worth zero.
        If the valuation date is before the swap payments start then this
        is the forward swap rate as it starts in the future. The swap rate
        is then a forward swap rate and so we use a forward discount
        factor. If the swap fixed leg has begun then we have a spot
        starting swap. The swap rate can also be calculated in a dual curve 
        approach but in this case the first fixing on the floating leg is
        needed. '''

        pv01 = self.pv01(valuationDate, discountCurve)

        if abs(pv01) < gSmall:
            raise FinError("PV01 is zero. Cannot compute swap rate.")

        if valuationDate < self._effectiveDate:
            df0 = discountCurve.df(self._effectiveDate)
        else:
            df0 = discountCurve.df(valuationDate)

        floatLegPV = 0.0

        if indexCurve is None:
            dfT = discountCurve.df(self._maturityDate)
            floatLegPV = (df0 - dfT) 
        else:
            floatLegPV = self._floatLeg.value(valuationDate,
                                              discountCurve,
                                              indexCurve, 
                                              firstFixing)

            floatLegPV /= self._fixedLeg._notional

        cpn = floatLegPV / pv01           
        return cpn
    
##########################################################################

    def cashSettledPV01(self,
                        valuationDate,
                        flatSwapRate,
                        frequencyType):
        ''' Calculate the forward value of an annuity of a forward starting
        swap using a single flat discount rate equal to the swap rate. This is
        used in the pricing of a cash-settled swaption in the FinIborSwaption
        class. This method does not affect the standard valuation methods.'''

        m = FinFrequency(frequencyType)

        if m == 0:
            raise FinError("Frequency cannot be zero.")

        ''' The swap may have started in the past but we can only value
        payments that have occurred after the valuation date. '''
        startIndex = 0
        while self._fixedLeg._paymentDates[startIndex] < valuationDate:
            startIndex += 1

        ''' If the swap has yet to settle then we do not include the
        start date of the swap as a coupon payment date. '''
        if valuationDate <= self._effectiveDate:
            startIndex = 1

        ''' Now PV fixed leg flows. '''
        flatPV01 = 0.0
        df = 1.0
        alpha = 1.0 / m

        for _ in self._fixedLeg._paymentDates[startIndex:]:
            df = df / (1.0 + alpha * flatSwapRate)
            flatPV01 += df * alpha

        return flatPV01

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
