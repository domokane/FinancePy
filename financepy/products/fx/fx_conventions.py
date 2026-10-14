##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...utils.date import Date
from dataclasses import dataclass

from ...utils.error import FinError
from ...utils.global_types import FXDeltaMethodTypes

########################################################################################


@dataclass(frozen=True)
class FXConvention:
    spot_lag: int
    settlement_calendars: tuple[str, ...]
    premium_currency: str | None = None
    delta_method: FXDeltaMethodTypes | None = None


G_FX_CONVENTIONS = {
    "EURUSD": FXConvention(
        spot_lag=2,
        settlement_calendars=("EUR", "USD"),
        premium_currency="USD",
        delta_method=FXDeltaMethodTypes.SPOT_DELTA,
    ),
    "USDJPY": FXConvention(
        spot_lag=2,
        settlement_calendars=("USD", "JPY"),
        premium_currency="USD",
        delta_method=FXDeltaMethodTypes.SPOT_DELTA_PREM_ADJ,
    ),
    "EURJPY": FXConvention(
        spot_lag=2,
        settlement_calendars=("EUR", "JPY"),
        premium_currency="EUR",
        delta_method=FXDeltaMethodTypes.SPOT_DELTA_PREM_ADJ,
    ),
    "USDCHF": FXConvention(
        spot_lag=2,
        settlement_calendars=("USD", "CHF"),
        premium_currency="USD",
        delta_method=FXDeltaMethodTypes.SPOT_DELTA_PREM_ADJ,
    ),
    "EURCHF": FXConvention(
        spot_lag=2,
        settlement_calendars=("EUR", "CHF"),
        premium_currency="EUR",
        delta_method=FXDeltaMethodTypes.SPOT_DELTA_PREM_ADJ,
    ),
    "GBPUSD": FXConvention(
        spot_lag=2,
        settlement_calendars=("GBP", "USD"),
        premium_currency="USD",
        delta_method=FXDeltaMethodTypes.SPOT_DELTA,
    ),
    "EURGBP": FXConvention(
        spot_lag=2,
        settlement_calendars=("EUR", "GBP"),
        premium_currency="EUR",
        delta_method=FXDeltaMethodTypes.SPOT_DELTA_PREM_ADJ,
    ),
    "AUDUSD": FXConvention(
        spot_lag=2,
        settlement_calendars=("AUD", "USD"),
        premium_currency="USD",
        delta_method=FXDeltaMethodTypes.SPOT_DELTA,
    ),
    "AUDJPY": FXConvention(
        spot_lag=2,
        settlement_calendars=("AUD", "JPY"),
        premium_currency="AUD",
        delta_method=FXDeltaMethodTypes.SPOT_DELTA_PREM_ADJ,
    ),
    "USDCAD": FXConvention(
        spot_lag=1,
        settlement_calendars=("USD", "CAD"),
        premium_currency="USD",
        delta_method=FXDeltaMethodTypes.SPOT_DELTA_PREM_ADJ,
    ),
    "USDBRL": FXConvention(
        spot_lag=2,
        settlement_calendars=("USD", "BRL"),
        premium_currency="USD",
        delta_method=FXDeltaMethodTypes.SPOT_DELTA_PREM_ADJ,
    ),
    "USDMXN": FXConvention(
        spot_lag=2,
        settlement_calendars=("USD", "MXN"),
        premium_currency="USD",
        delta_method=FXDeltaMethodTypes.SPOT_DELTA_PREM_ADJ,
    ),
}

# Additional Cboe FX spot pairs.
# Existing FX_CONVENTIONS entries are retained.
# Option conventions are deliberately left as None.

ADDITIONAL_FX_SPOT_LAGS = {
    # AUD crosses
    "AUDCAD": 2,
    "AUDCHF": 2,
    "AUDHKD": 2,
    "AUDNZD": 2,

    # CAD crosses
    "CADCHF": 2,
    "CADJPY": 2,
    "CADMXN": 2,

    # CHF crosses
    "CHFJPY": 2,
    "CHFNOK": 2,
    "CHFSEK": 2,

    # EUR crosses
    "EURAUD": 2,
    "EURCAD": 2,
    "EURCNH": 2,
    "EURCZK": 2,
    "EURDKK": 2,
    "EURHKD": 2,
    "EURHUF": 2,
    "EURMXN": 2,
    "EURNOK": 2,
    "EURNZD": 2,
    "EURPLN": 2,
    "EURRON": 2,
    "EURSEK": 2,
    "EURSGD": 2,
    "EURTRY": 2,
    "EURZAR": 2,

    # GBP crosses
    "GBPAUD": 2,
    "GBPCAD": 2,
    "GBPCHF": 2,
    "GBPCZK": 2,
    "GBPHUF": 2,
    "GBPJPY": 2,
    "GBPMXN": 2,
    "GBPNOK": 2,
    "GBPNZD": 2,
    "GBPPLN": 2,
    "GBPSEK": 2,

    # Other crosses
    "HKDJPY": 2,
    "NOKSEK": 2,
    "NZDJPY": 2,
    "NZDUSD": 2,
    "ZARJPY": 2,

    # USD crosses
    "USDAED": 2,
    "USDBHD": 2,
    "USDCNH": 2,
    "USDCZK": 2,
    "USDDKK": 2,
    "USDHKD": 2,
    "USDHUF": 2,
    "USDILS": 2,
    "USDJOD": 2,
    "USDKWD": 2,
    "USDMAD": 2,
    "USDNOK": 2,
    "USDOMR": 2,
    "USDPLN": 2,
    "USDQAR": 2,
    "USDRON": 2,
    "USDSAR": 2,
    "USDSEK": 2,
    "USDSGD": 2,
    "USDTHB": 2,
    "USDTND": 2,
    "USDTRY": 1,
    "USDZAR": 2,
}

G_FX_CONVENTIONS.update({
    pair: FXConvention(
        spot_lag=lag,
        settlement_calendars=(pair[:3], pair[3:]),
    )
    for pair, lag in ADDITIONAL_FX_SPOT_LAGS.items()
})


def get_fx_convention(currency_pair: str) -> FXConvention:
    """Return the market conventions for an FX currency pair."""

    pair = currency_pair.upper()

    if pair not in G_FX_CONVENTIONS:
        raise FinError(f"Unsupported FX currency pair: {pair}")

    return G_FX_CONVENTIONS[pair]


def get_fx_spot_lag(currency_pair: str) -> int:
    """Return the conventional FX spot settlement lag."""
    return get_fx_convention(currency_pair).spot_lag


prem_currency = {
    pair: convention.premium_currency
    for pair, convention in G_FX_CONVENTIONS.items()
    if convention.premium_currency is not None
}

delta_convention = {
    pair: convention.delta_method
    for pair, convention in G_FX_CONVENTIONS.items()
    if convention.delta_method is not None
}


if 1 == 1:
    from ...utils.currency import g_ccy_names

    for pair, convention in G_FX_CONVENTIONS.items():

        assert len(pair) == 6

        ccy1 = pair[:3]
        ccy2 = pair[3:]

        assert ccy1 in g_ccy_names
        assert ccy2 in g_ccy_names
        assert ccy1 != ccy2

        assert convention.spot_lag >= 0

        assert convention.settlement_calendars == (ccy1, ccy2)

        if convention.premium_currency is not None:
            assert convention.premium_currency in (ccy1, ccy2)
