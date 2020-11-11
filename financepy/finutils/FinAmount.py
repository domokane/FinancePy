##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from .FinHelperFunctions import checkArgumentTypes
from .FinCurrency import FinCurrencyTypes
from .FinMath import ONE_MILLION

###############################################################################


class FinAmount(object):
    ''' A FinAmount is a holder for an amount in a specific currency. '''

    def __init__(self,
                 notional: float = ONE_MILLION,
                 currencyType: FinCurrencyTypes = FinCurrencyTypes.NONE):
        ''' Create FinAmount object. '''

        checkArgumentTypes(self.__init__, locals())

        self._notional = notional
        self._currencyType = currencyType


    def __repr__(self):
        ''' Print out the details of the schedule and the actual dates. This
        can be used for providing transparency on schedule calculations. '''

        s = ""
        if self._currencyType != FinCurrencyTypes.NONE:
            s += self._currencyType.name
            s += " "

        s += '{:,.2f}'.format(self._notional)
        
        return s

    def amount(self):
        return self._notional

    def _print(self):
        ''' Print out the details of the schedule and the actual dates. This
        can be used for providing transparency on schedule calculations. '''
        print(self)

###############################################################################


