##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO: Add functionality around settlement
# TODO: Write test function
# TODO: Handle 1 month futures contracts

import numpy as np

from ...finutils.FinError import FinError
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinMath import ONE_MILLION
from ...finutils.FinDate import FinDate

from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...products.rates.FinIborFRA import FinIborFRA

###############################################################################


class FinIborFuture(object):
    ''' Class for managing short term interest rate futures contracts. '''

    # Reference
    # https://www.cmegroup.com/education/files/eurodollar-futures-the-basics-file01.pdf

    def __init__(self,
                 todayDate: FinDate,
                 futureNumber: int,  # The number of the future after todayDate
                 futureTenor: str = "3M",  # '1M', '2M', '3M'
                 accrualType: FinDayCountTypes = FinDayCountTypes.ACT_360,
                 contractSize: float = ONE_MILLION):
        ''' Create an interest rate futures contract which has the same
        conventions as those traded on the CME. The current date, the tenor of
        the future, the number of the future and the accrual convention and
        the contract size should be provided. '''

        checkArgumentTypes(self.__init__, locals())

        if futureNumber < 1:
            raise FinError("Future number must be 1 or more")

        if futureTenor != "3M" and futureTenor != "3m":
            raise FinError("Only 3M IMM futures handled currently.")

        self._deliveryDate = todayDate.nextIMMDate()

        for iFut in range(0, futureNumber - 1):
            self._deliveryDate = self._deliveryDate.nextIMMDate()

        self._endOfInterestPeriod = self._deliveryDate.nextIMMDate()

        self._lastTradingDate = self._deliveryDate.addDays(-2)
        self._accrualType = accrualType
        self._contractSize = contractSize

###############################################################################

    def toFRA(self, futuresPrice, convexity):
        ''' Convert the futures contract to a FinIborFRA object so it can be
        used to boostrap a Ibor curve. For this we need to adjust the futures
        rate using the convexity correction. '''

        fraRate = self.FRARate(futuresPrice, convexity)

        fra = FinIborFRA(self._deliveryDate,
                          self._endOfInterestPeriod,
                          fraRate,
                          self._accrualType,
                          notional=self._contractSize,
                          payFixedRate=False)

        return fra

###############################################################################

    def futuresRate(self, futuresPrice):
        ''' Calculate implied futures rate from the futures price.'''
        futuresRate = (100.0 - futuresPrice) / 100.0
        return futuresRate

###############################################################################

    def FRARate(self, futuresPrice, convexity):
        ''' Convert futures price and convexity to a FRA rate using the BBG
        negative convexity (in percent). This is then divided by 100 before
        being added to the futures rate. '''

        futRate = (100.0 - futuresPrice) / 100.0

        if convexity < 0:
            fraRate = futRate + convexity/100.0
        else:
            fraRate = futRate - convexity/100.0

        return fraRate

###############################################################################

    def convexity(self, valuationDate, volatility, meanReversion):
        ''' Calculation of the convexity adjustment between FRAs and interest
        rate futures using the Hull-White model as described in technical note
        in link below:
        http://www-2.rotman.utoronto.ca/~hull/TechnicalNotes/TechnicalNote1.pdf

        NOTE THIS DOES NOT APPEAR TO AGREE WITH BLOOMBERG!! INVESTIGATE.
        '''

        a = meanReversion
        t0 = 0.0
        t1 = (self._lastTradingDate - valuationDate) / gDaysInYear
        t2 = (self._endOfInterestPeriod - valuationDate) / gDaysInYear

        # Hull White model for short rate dr = (theta(t)-ar) dt + sigma * dz
        # This reduces to Ho-Lee when a = 0 so to avoid divergences I provide
        # this numnerical limit
        if abs(a) > 1e-10:

            bt1t2 = (1.0 - np.exp(-a * (t2 - t1))) / a
            bt0t1 = (1.0 - np.exp(-a * (t1 - t0))) / a
            w = 1.0 - np.exp(-2.0 * a * t1)
            term = bt1t2 * w + 2.0 * a * (bt0t1**2)
            c = bt1t2 * (volatility**2) * term / (t2 - t1) / 4.0 / a

        else:
            c = t1 * t2 * (volatility**2) / 2.0

        return c

##########################################################################

    def __repr__(self):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("LAST TRADING DATE", self._lastTradingDate)
        s += labelToString("DELIVERY DATE", self._deliveryDate)
        s += labelToString("END INTEREST PERIOD", self._endOfInterestPeriod)
        s += labelToString("ACCRUAL TYPE", self._accrualType)
        s += labelToString("CONTRACT SIZE", self._contractSize)
        return s

##########################################################################
