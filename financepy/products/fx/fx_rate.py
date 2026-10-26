##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...utils.fx_conventions import G_FX_CONVENTIONS
from ...utils.calendar import Calendar
from ...utils.date import Date
from ...utils.currency import get_currency_calendar_type
from ...utils.check_values import check_currency_pair
from ...utils.error import FinError
import math

# DO NOT USE


# class FXRate:
#     """Exchange rate quoted in units of the second currency
#     per unit of the first currency.

#     Can represent a spot rate, forward rate, or contractual
#     forward strike. The currency-pair order determines the
#     quotation direction.

#     For example, FXRate("EURUSD", 1.1206) represents
#     1 EUR = 1.1206 USD.
#     """

#     def __init__(self, currency_pair: str, rate: float) -> None:

#         currency_pair = currency_pair.upper()

#         check_currency_pair(currency_pair)

#         if not math.isfinite(rate) or rate <= 0.0:
#             raise FinError("FX rate must be positive and finite.")

#         self.currency_pair = currency_pair
#         self.ccy1 = currency_pair[:3]
#         self.ccy2 = currency_pair[3:]
#         self.value = rate

#     def currency_pair_names(self):
#         """Return the currency names."""

#         return [self.ccy1, self.ccy2]

#     def inverse(self) -> "FXRate":
#         """Return the reciprocal FX quotation."""

#         inverse_pair = self.ccy2 + self.ccy1

#         return FXRate(inverse_pair, 1.0 / self.rate)

#     def bumnp(self, bump_size) -> "FXRate":
#         """Return the reciprocal FX quotation."""

#         bumped_rate = self.value + bump_size

#         return FXRate(self.currency_pair, bumped_rate)

#     def spot_date(self, trade_dt: Date) -> Date:
#         """Return the conventional FX spot settlement date."""

#         cal_type_1 = get_currency_calendar_type(self.ccy1)
#         cal_type_2 = get_currency_calendar_type(self.ccy2)

#         calendar = Calendar([cal_type_1, cal_type_2])

#         spot_lag = G_FX_CONVENTIONS[self.currency_pair].spot_lag

#         spot_dt = calendar.add_business_days(trade_dt, spot_lag)

#         return spot_dt

#     def __sub__(self, other: "FXRate") -> float:
#         if not isinstance(other, FXRate):
#             return NotImplemented

#         if self.currency_pair != other.currency_pair:
#             raise FinError("Cannot subtract FX rates with different currency pairs.")

#         return FXRate(self.currency_pair, self.value - other.value)

#     def __repr__(self) -> str:
#         return (
#             f"FXRate(currency_pair='{self.currency_pair}', "
#             f"rate={self.value})"
#         )
