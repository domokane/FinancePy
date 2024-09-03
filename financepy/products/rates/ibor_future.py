##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO: Add functionality around settlement
# TODO: Write test function
# TODO: Handle 1 month futures contracts

import numpy as np

from ...utils.error import FinError
from ...utils.day_count import DayCountTypes
from ...utils.global_vars import g_days_in_year
from ...utils.math import ONE_MILLION
from ...utils.date import Date

from ...utils.helpers import label_to_string, check_argument_types
from ...products.rates.ibor_fra import IborFRA

###############################################################################


class IborFuture:
    """Class for managing short term interest rate futures contracts."""

    # Reference
    # https://www.cmegroup.com/education/files/eurodollar-futures-the-basics-file01.pdf

    def __init__(
        self,
        today_dt: Date,
        future_number: int,  # The number of the future after today_dt
        futureTenor: str = "3M",  # '1M', '2M', '3M'
        dc_type: DayCountTypes = DayCountTypes.ACT_360,
        contract_size: float = ONE_MILLION,
    ):
        """Create an interest rate futures contract which has the same
        conventions as those traded on the CME. The current _dt, the tenor of
        the future, the number of the future and the accrual convention and
        the contract size should be provided."""

        check_argument_types(self.__init__, locals())

        if future_number < 1:
            raise FinError("Future number must be 1 or more")

        if futureTenor != "3M" and futureTenor != "3m":
            raise FinError("Only 3M IMM futures handled currently.")

        self.delivery_dt = today_dt.next_imm_date()

        for iFut in range(0, future_number - 1):
            self.delivery_dt = self.delivery_dt.next_imm_date()

        self.end_of_interest_period = self.delivery_dt.next_imm_date()

        self.last_trading_dt = self.delivery_dt.add_days(-2)
        self.dc_type = dc_type
        self.contract_size = contract_size

    ###########################################################################

    def to_fra(self, futures_price, convexity):
        """Convert the futures contract to a IborFRA object so it can be
        used to boostrap a Ibor curve. For this we need to adjust the futures
        rate using the convexity correction."""

        fra_rate = self.fra_rate(futures_price, convexity)

        fra = IborFRA(
            self.delivery_dt,
            self.end_of_interest_period,
            fra_rate,
            self.dc_type,
            notional=self.contract_size,
            pay_fixed_rate=False,
        )

        return fra

    ###########################################################################

    def futures_rate(self, futures_price):
        """Calculate implied futures rate from the futures price."""
        futures_rate = (100.0 - futures_price) / 100.0
        return futures_rate

    ###########################################################################

    def fra_rate(self, futures_price, convexity):
        """Convert futures price and convexity to a FRA rate using the BBG
        negative convexity (in percent). This is then divided by 100 before
        being added to the futures rate."""

        futRate = (100.0 - futures_price) / 100.0

        if convexity < 0:
            fra_rate = futRate + convexity / 100.0
        else:
            fra_rate = futRate - convexity / 100.0

        return fra_rate

    ###########################################################################

    def convexity(self, value__dt, volatility, mean_reversion):
        """Calculation of the convexity adjustment between FRAs and interest
        rate futures using the Hull-White model as described in technical note
        in link below:
        http://www-2.rotman.utoronto.ca/~hull/TechnicalNotes/TechnicalNote1.pdf

        NOTE THIS DOES NOT APPEAR TO AGREE WITH BLOOMBERG!! INVESTIGATE.
        """

        a = mean_reversion
        t0 = 0.0
        t1 = (self.last_trading_dt - value__dt) / g_days_in_year
        t2 = (self.end_of_interest_period - value__dt) / g_days_in_year

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

    ###########################################################################

    def __repr__(self):
        """Print a list of the unadjusted coupon payment _dts used in
        analytic calculations for the bond."""
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("LAST TRADING DATE", self.last_trading_dt)
        s += label_to_string("DELIVERY DATE", self.delivery_dt)
        s += label_to_string(
            "END INTEREST PERIOD", self.end_of_interest_period
        )
        s += label_to_string("DAY COUNT TYPE", self.dc_type)
        s += label_to_string("CONTRACT SIZE", self.contract_size)
        return s


###############################################################################
