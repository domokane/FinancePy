##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinCalendar import FinCalendarTypes, FinDateGenRuleTypes
from ...finutils.FinCalendar import FinCalendar, FinBusDayAdjustTypes
from ...finutils.FinHelperFunctions import checkArgumentTypes, labelToString
from ...finutils.FinMath import ONE_MILLION
from ...finutils.FinGlobalTypes import FinSwapTypes
from ...market.curves.FinDiscountCurve import FinDiscountCurve

from .FinFloatLeg import FinFloatLeg

###############################################################################


class FinIborOIS(object):
    ''' Class for managing an Ibor-OIS basis swap contract. This is a
    contract in which a floating leg with one LIBOR tenor is exchanged for a 
    floating leg payment of an overnight index swap. There is no exchange of
    par. The contract is entered into at zero initial cost. The contract lasts
    from a start date to a specified maturity date.
    
    The value of the contract is the NPV of the two coupon streams. Discounting
    is done on a supplied discount curve which is separate from the curves from
    which the implied index rates are extracted. '''
    
    def __init__(self,
                 effectiveDate: FinDate,  # Date interest starts to accrue
                 terminationDateOrTenor: (FinDate, str),  # Date contract ends
                 iborType: FinSwapTypes,
                 iborFreqType: FinFrequencyTypes = FinFrequencyTypes.QUARTERLY,
                 iborDayCountType: FinDayCountTypes  = FinDayCountTypes.THIRTY_E_360,
                 iborSpread: float = 0.0,
                 oisFreqType: FinFrequencyTypes = FinFrequencyTypes.QUARTERLY,
                 oisDayCountType: FinDayCountTypes = FinDayCountTypes.THIRTY_E_360,
                 oisSpread: float = 0.0,
                 oisPaymentLag: int = 0,
                 notional: float = ONE_MILLION,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD):
        ''' Create a Ibor basis swap contract giving the contract start
        date, its maturity, frequency and day counts on the two floating 
        legs and notional. The floating leg parameters have default
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

        oisType = FinSwapTypes.PAY
        if iborType == FinSwapTypes.PAY:
            oisType = FinSwapTypes.RECEIVE
        
        principal = 0.0

        self._floatIborLeg = FinFloatLeg(effectiveDate,
                                         self._terminationDate,
                                         iborType,
                                         iborSpread,
                                         iborFreqType,
                                         iborDayCountType,
                                         notional,
                                         principal,
                                         0,
                                         calendarType,
                                         busDayAdjustType,
                                         dateGenRuleType)

        self._floatOISLeg = FinFloatLeg(effectiveDate,
                                        self._terminationDate,
                                        oisType,
                                        oisSpread,
                                        oisFreqType,
                                        oisDayCountType,
                                        notional,
                                        principal,
                                        oisPaymentLag,
                                        calendarType,
                                        busDayAdjustType,
                                        dateGenRuleType)

###############################################################################

    def value(self,
              valuationDate: FinDate,
              discountCurve: FinDiscountCurve,
              indexIborCurve: FinDiscountCurve = None,
              indexOISCurve: FinDiscountCurve = None,
              firstFixingRateLeg1=None,
              firstFixingRateLeg2=None):
        ''' Value the interest rate swap on a value date given a single Ibor
        discount curve and an index curve for the Ibors on each swap leg. '''

        if indexIborCurve is None:
            indexIborCurve = discountCurve

        if indexOISCurve is None:
            indexOISCurve = discountCurve

        floatIborLegValue = self._floatIborLeg.value(valuationDate,
                                                     discountCurve, 
                                                     indexIborCurve, 
                                                     firstFixingRateLeg1)

        floatOISLegValue = self._floatOISLeg.value(valuationDate,
                                                   discountCurve,
                                                   indexOISCurve,
                                                   firstFixingRateLeg2)

        value = floatIborLegValue + floatOISLegValue
        return value

###############################################################################

    def printFlows(self):
        ''' Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows. '''

        self._floatIborLeg.printPayments()
        self._floatOISLeg.printPayments()

##########################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += self._floatIborLeg.__repr__()
        s += "\n"
        s += self._floatOISLeg.__repr__()
        return s

###############################################################################

    def _print(self):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''
        print(self)

###############################################################################
