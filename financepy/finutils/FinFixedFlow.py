##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from .FinError import FinError
from .FinDate import FinDate
from .FinHelperFunctions import labelToString
from .FinHelperFunctions import checkArgumentTypes
from .FinDayCount import FinDayCount, FinDayCountTypes

###############################################################################

# UNFINISHED

class FinIborFlow(object):
    ''' A FinIborFlow is a quantity of cash paid at some time where the
    quantity is based on an Ibor interest rate, an accrual period, a notional
    size and a currency. '''

    def __init__(self,
                 startAccruedDate: FinDate,
                 endAccruedDate: FinDate,
                 paymentDate: FinDate,
                 dayCountType: FinDayCountTypes,
                 notional: float):
        ''' Create FinFixedFlow object. '''

        checkArgumentTypes(self.__init__, locals())

        if startAccruedDate >= endAccruedDate:
            raise FinError("Start accrued date must precede end accrued date.")

        # validation complete
        self._startAccruedDate = startAccruedDate
        self._endAccruedDate = endAccruedDate
        self._paymentDate = paymentDate
        self._notional = notional
        self._dayCountType = dayCountType


    def amount(self, indexCurve):
        
        rate = indexCurve.df()
        dayCounter = FinDayCount(self._dayCountType)
        self._yearFrac = dayCounter.yearFrac(self._startAccruedDate,
                                             self._endAccruedDate)
        self._cashflow = self._yearFrac * self._rate * self._notional
        return self._cashflow

##############################################################################

    def __repr__(self):
        ''' Print out the details of the schedule and the actual dates. This
        can be used for providing transparency on schedule calculations. '''

        s = labelToString("START ACCRUED DATE", self._startAccruedDate)
        s += labelToString("END ACCRUED DATE", self._endAccruedDate)
        s += labelToString("PAYMENT DATE", self._paymentDate)
        s += labelToString("RATE", self._rate)
        s += labelToString("NOTIONAL", self._notional)
        s += labelToString("DAY COUNT TYPE", self._dayCountType)
        return s

###############################################################################

    def _print(self):
        ''' Print out the details of the schedule and the actual dates. This
        can be used for providing transparency on schedule calculations. '''
        print(self)

###############################################################################

