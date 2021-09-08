##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO: Add functionality around settlement
# TODO: Write test function
# TODO: Handle 1 month futures contracts

import numpy as np

from ...utils.error import FinError
from ...utils.day_count import DayCountTypes
from ...utils.global_vars import gDaysInYear
from ...utils.math import ONE_MILLION
from ...utils.date import Date

from ...utils.helpers import label_to_string, check_argument_types
from ...products.rates.ibor_fra import IborFRA

###############################################################################


class IborFuture:
    """ Class for managing short term interest rate futures contracts. """

    # Reference
    # https://www.cmegroup.com/education/files/eurodollar-futures-the-basics-file01.pdf

    def __init__(self,
                 todayDate: Date,
                 futureNumber: int,  # The number of the future after todayDate
                 futureTenor: str = "3M",  # '1M', '2M', '3M'
                 accrual_type: DayCountTypes = DayCountTypes.ACT_360,
                 contract_size: float = ONE_MILLION):
        """ Create an interest rate futures contract which has the same
        conventions as those traded on the CME. The current date, the tenor of
        the future, the number of the future and the accrual convention and
        the contract size should be provided. """

        check_argument_types(self.__init__, locals())

        if futureNumber < 1:
            raise FinError("Future number must be 1 or more")

        if futureTenor != "3M" and futureTenor != "3m":
            raise FinError("Only 3M IMM futures handled currently.")

        self._delivery_date = todayDate.next_imm_date()

        for iFut in range(0, futureNumber - 1):
            self._delivery_date = self._delivery_date.next_imm_date()

        self._endOfInterestPeriod = self._delivery_date.next_imm_date()

        self._lastTradingDate = self._delivery_date.add_days(-2)
        self._accrual_type = accrual_type
        self._contract_size = contract_size

###############################################################################

    def to_fra(self, futures_price, convexity):
        """ Convert the futures contract to a IborFRA object so it can be
        used to boostrap a Ibor curve. For this we need to adjust the futures
        rate using the convexity correction. """

        fraRate = self.fra_rate(futures_price, convexity)

        fra = IborFRA(self._delivery_date,
                      self._endOfInterestPeriod,
                      fraRate,
                      self._accrual_type,
                      notional=self._contract_size,
                      payFixedRate=False)

        return fra

###############################################################################

    def futures_rate(self, futures_price):
        """ Calculate implied futures rate from the futures price."""
        futuresRate = (100.0 - futures_price) / 100.0
        return futuresRate

###############################################################################

    def fra_rate(self, futures_price, convexity):
        """ Convert futures price and convexity to a FRA rate using the BBG
        negative convexity (in percent). This is then divided by 100 before
        being added to the futures rate. """

        futRate = (100.0 - futures_price) / 100.0

        if convexity < 0:
            fraRate = futRate + convexity/100.0
        else:
            fraRate = futRate - convexity/100.0

        return fraRate

###############################################################################

    def convexity(self, valuation_date, volatility, mean_reversion):
        """ Calculation of the convexity adjustment between FRAs and interest
        rate futures using the Hull-White model as described in technical note
        in link below:
        http://www-2.rotman.utoronto.ca/~hull/TechnicalNotes/TechnicalNote1.pdf

        NOTE THIS DOES NOT APPEAR TO AGREE WITH BLOOMBERG!! INVESTIGATE.
        """

        a = mean_reversion
        t0 = 0.0
        t1 = (self._lastTradingDate - valuation_date) / gDaysInYear
        t2 = (self._endOfInterestPeriod - valuation_date) / gDaysInYear

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
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("LAST TRADING DATE", self._lastTradingDate)
        s += label_to_string("DELIVERY DATE", self._delivery_date)
        s += label_to_string("END INTEREST PERIOD", self._endOfInterestPeriod)
        s += label_to_string("ACCRUAL TYPE", self._accrual_type)
        s += label_to_string("CONTRACT SIZE", self._contract_size)
        return s

##########################################################################
