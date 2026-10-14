from .fx_forward import FXForward
from .fx_conventions import get_fx_spot_lag
from .fx_rate import FXRate
from ...utils.currency import ccy_names
from ...utils.date import Date


class FXSpot(FXForward):
    """Deliverable FX spot transaction.

    The holder receives the foreign currency and pays the
    domestic currency at the agreed exchange rate.

    The settlement date is determined from the trade date
    using the currency-pair spot settlement conventions.
    """

    def __init__(
        self,
        trade_dt: Date,
        strike_fx_rate: FXRate,
        notional: float,
        notional_currency: str,
        buy_currency: str
    ) -> None:

        settlement_dt = self.spot_date(
            trade_dt,
            strike_fx_rate.currency_pair,
        )

        super().__init__(
            settlement_dt=settlement_dt,
            strike_fx_rate=strike_fx_rate,
            notional=notional,
            notional_currency=notional_currency,
            buy_currency=buy_currency
        )

        self.trade_dt = trade_dt

    def spot_date(self, trade_dt: Date) -> Date:
        """Return the conventional FX spot settlement date.

        Uses weekend-only business-day adjustment.
        Currency holidays are not yet incorporated.
        """

        spot_lag = get_fx_spot_lag(self.currency_pair)

        ccy_pair_names = self.strike_fx_rate.ccy_pair_names()

        asjusted_dt = trade
        return trade_dt.add_weekdays(spot_lag)
