##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...utils.error import FinError
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.global_vars import G_DAYS_IN_YEAR
from ...utils.math import ONE_MILLION
from ...utils.date import Date

from ...utils.helpers import label_to_string, check_argument_types
from ...products.rates.ibor_fra import IborFRA

########################################################################################


class IborFuture:
    """Class for managing short term interest rate futures contracts."""

    # Reference
    # https://www.cmegroup.com/education/files/eurodollar-futures-the-basics-file01.pdf

    def __init__(
        self,
        today_dt: Date,
        future_number: int,  # The number of the future after today_dt
        future_tenor: str = "3M",  # '1M', '3M'
        accrual_dc_type: DayCountTypes = DayCountTypes.ACT_360,
        contract_size: float = ONE_MILLION,
    ):
        """Create an interest rate futures contract which has the same
        conventions as those traded on the CME. The current _dt, the tenor of
        the future, the number of the future and the accrual convention and
        the contract size should be provided."""

        check_argument_types(self.__init__, locals())

        if future_number < 1:
            raise FinError("Future number must be 1 or more")

        future_tenor = future_tenor.upper()

        if future_tenor not in ["1M", "3M"]:
            raise FinError("Only 1M and 3M IMM futures handled currently.")

        self.future_tenor = future_tenor
        self.delivery_dt = self._first_delivery_dt(today_dt, future_tenor)

        for _ in range(0, future_number - 1):
            self.delivery_dt = self._next_delivery_dt(
                self.delivery_dt, future_tenor
            )

        self.end_of_interest_period = self._next_delivery_dt(
            self.delivery_dt, future_tenor
        )
        self.last_trading_dt = self.delivery_dt.add_days(-2)
        self.accrual_dc_type = accrual_dc_type
        self.contract_size = contract_size

    ###########################################################################

    @staticmethod
    def _next_monthly_imm_date(dt: Date):
        """Return the next third Wednesday after dt, for any month."""

        m_imm = dt.m
        y_imm = dt.y
        d_imm = dt.third_wednesday_of_month(m_imm, y_imm)

        if dt.d >= d_imm:
            next_month_dt = dt.add_months(1)
            m_imm = next_month_dt.m
            y_imm = next_month_dt.y
            d_imm = next_month_dt.third_wednesday_of_month(m_imm, y_imm)

        return Date(d_imm, m_imm, y_imm)

    ###########################################################################

    @staticmethod
    def _first_delivery_dt(today_dt: Date, future_tenor: str):
        """Return the first delivery date after today for the tenor."""

        if future_tenor == "1M":
            return IborFuture._next_monthly_imm_date(today_dt)

        return today_dt.next_imm_date()

    ###########################################################################

    @staticmethod
    def _next_delivery_dt(delivery_dt: Date, future_tenor: str):
        """Return the next delivery date after an existing delivery date."""

        if future_tenor == "1M":
            return IborFuture._next_monthly_imm_date(delivery_dt)

        return delivery_dt.next_imm_date()

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
            self.accrual_dc_type,
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

    def accrual_factor(self):
        """Return the accrual factor for the futures interest period."""

        dc = DayCount(self.accrual_dc_type)
        return dc.year_frac(
            self.delivery_dt, self.end_of_interest_period
        )[0]

    ###########################################################################

    def basis_point_value(self):
        """Return the cash value of one basis point move in futures price."""

        return self.contract_size * self.accrual_factor() / 10000.0

    ###########################################################################

    def settlement_amount(
        self,
        futures_price: float,
        settlement_price: float,
        num_contracts: float = 1.0,
    ):
        """Return the cash settlement amount for a futures price move.

        The amount is for a long position. Use a negative number of contracts
        for a short position. Prices are quoted in futures price points, for
        example 97.50.
        """

        price_change = settlement_price - futures_price
        return (
            num_contracts
            * self.contract_size
            * self.accrual_factor()
            * price_change
            / 100.0
        )

    ###########################################################################

    def fra_rate(self, futures_price, convexity):
        """Convert futures price and convexity to a FRA rate using the BBG
        negative convexity (in percent). This is then divided by 100 before
        being added to the futures rate."""

        futures_rate = (100.0 - futures_price) / 100.0

        if convexity < 0:
            fra_rate = futures_rate + convexity / 100.0
        else:
            fra_rate = futures_rate - convexity / 100.0

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
        t1 = (self.last_trading_dt - value__dt) / G_DAYS_IN_YEAR
        t2 = (self.end_of_interest_period - value__dt) / G_DAYS_IN_YEAR

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
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("LAST TRADING DATE", self.last_trading_dt)
        s += label_to_string("DELIVERY DATE", self.delivery_dt)
        s += label_to_string(
            "END INTEREST PERIOD", self.end_of_interest_period
        )
        s += label_to_string("DC_TYPE", self.accrual_dc_type)
        s += label_to_string("CONTRACT SIZE", self.contract_size)
        return s


########################################################################################
