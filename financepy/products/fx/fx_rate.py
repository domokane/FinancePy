##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import math

from ...utils.error import FinError
from ...utils.check_values import check_currency_pair
from ...utils.currency import get_currency_calendar_type
from ...utils.date import Date
from ...utils.calendar import Calendar
from .fx_conventions import G_FX_CONVENTIONS


class FXRate:
    """Exchange rate quoted in units of the second currency
    per unit of the first currency.

    Can represent a spot rate, forward rate, or contractual
    forward strike. The currency-pair order determines the
    quotation direction.

    For example, FXRate("EURUSD", 1.1206) represents
    1 EUR = 1.1206 USD.
    """

    def __init__(self, ccy_pair: str, rate: float) -> None:

        ccy_pair = ccy_pair.upper()

        check_currency_pair(ccy_pair)

        if not math.isfinite(rate) or rate <= 0.0:
            raise FinError("FX rate must be positive and finite.")

        self.ccy_pair = ccy_pair
        self.ccy1 = ccy_pair[:3]
        self.ccy2 = ccy_pair[3:]
        self.rate = rate

    def ccy_pair_names(self):
        """Return the currency names."""

        return [self.ccy1, self.ccy2]

    def inverse(self) -> "FXRate":
        """Return the reciprocal FX quotation."""

        inverse_pair = self.ccy2 + self.ccy1

        return FXRate(inverse_pair, 1.0 / self.rate)

    def spot_date(self, trade_dt: Date) -> Date:
        """Return the conventional FX spot settlement date."""

        cal_type_1 = get_currency_calendar_type(self.ccy1)
        cal_type_2 = get_currency_calendar_type(self.ccy2)

        calendar = Calendar([cal_type_1, cal_type_2])

        spot_lag = G_FX_CONVENTIONS[self.ccy_pair].spot_lag

        spot_dt = calendar.add_business_days(trade_dt, spot_lag)

        return spot_dt

    def __sub__(self, other: "FXRate") -> float:
        if not isinstance(other, FXRate):
            return NotImplemented

        if self.ccy_pair != other.ccy_pair:
            raise FinError("Cannot subtract FX rates with different currency pairs.")

        return FXRate(self.ccy_pair, self.rate - other.rate)

    def __repr__(self) -> str:
        return (
            f"FXRate(currency_pair='{self.ccy_pair}', "
            f"rate={self.rate})"
        )
