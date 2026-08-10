##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from ...utils.error import FinError
from ...utils.global_types import FXDeltaMethodTypes

# Non exhaustive list of country codes and currency names

ccyNames = {
    "AED": ("UNITED ARAB EMIRATES", "DIRHAM"),
    "AUD": ("AUSTRALIA", "DOLLAR"),
    "BRL": ("BRAZIL", "REAL"),
    "CAD": ("CANADA", "DOLLAR"),
    "CHF": ("SWITZERLAND", "FRANC"),
    "CLP": ("CHILE", "PESO"),
    "CNY": ("CHINA", "RENMIMBI"),
    "COP": ("COLUMBIA", "PESO"),
    "DKK": ("DENMARK", "KRONE"),
    "EUR": ("EUROZONE", "EURO"),
    "GBP": ("UK", "POUND"),
    "HKD": ("HONG KONG", "DOLLAR"),
    "HUF": ("HUNGARY", "FORINT"),
    "IDR": ("INDONESIA", "RUPIAH"),
    "INR": ("INDIA", "RUPEE"),
    "ILS": ("ISRAEL", "SHEKEL"),
    "JPY": ("JAPAN", "YEN"),
    "KRW": ("KOREAN", "WON"),
    "MYR": ("MALYSIA", "RINGIT"),
    "MXN": ("MEXICO", "PESO"),
    "NOK": ("NORWAY", "KRONE"),
    "NZD": ("NEW ZEALAND", "DOLLAR"),
    "PHP": ("PHILIPPINES", "PESO"),
    "PLN": ("POLAND", "ZLOTY"),
    "RON": ("ROMANIA", "LEU"),
    "RUB": ("RUSSIA", "RUBLE"),
    "SAR": ("SAUDI ARABIA", "RIYAL"),
    "SEK": ("SWEDEN", "KRONA"),
    "SGD": ("SINGAPORE", "DOLLAR"),
    "THB": ("THAILAND", "BHAT"),
    "TRY": ("TURKEY", "LIRA"),
    "TWD": ("TAIWAN", "DOLLAR"),
    "USD": ("US", "DOLLAR"),
    "ZAR": ("SOUTH AFRICA", "RAND"),
}

ccyQuotes = {
    "EURUSD": (1, 1, 1),
    "USDJPY": (1, 1, 1),
    "EURJPY": (1, 1, 1),
    "GBPUSD": (1, 1, 1),
    "EURGBP": (1, 1, 1),
    "USDCHF": (1, 1, 1),
    "AUDUSD": (1, 1, 1),
    "NZDUSD": (1, 1, 1),
    "USDCAD": (1, 1, 1),
    "EURNOK": (1, 1, 1),
    "EURSEK": (1, 1, 1),
    "EURDKK": (1, 1, 1),
    "EURHUF": (1, 1, 1),
    "EURPLN": (1, 1, 1),
    "USDTRY": (1, 1, 1),
    "USDZAR": (1, 1, 1),
    "USDMXN": (1, 1, 1),
    "USDBRL": (1, 1, 1),
    "USDSGD": (1, 1, 1),
}

########################################################################################


########################################################################################


prem_currency = {
    "EURUSD": "USD",
    "USDJPY": "USD",
    "EURJPY": "EUR",
    "USDCHF": "USD",
    "EURCHF": "EUR",
    "GBPUSD": "USD",
    "EURGBP": "EUR",
    "AUDUSD": "USD",
    "AUDJPY": "AUD",
    "USDCAD": "USD",
    "USDBRL": "USD",
    "USDMXN": "USD",
}

delta_convention = {
    "EURUSD": FXDeltaMethodTypes.SPOT_DELTA,
    "USDJPY": FXDeltaMethodTypes.SPOT_DELTA_PREM_ADJ,
    "EURJPY": FXDeltaMethodTypes.SPOT_DELTA_PREM_ADJ,
    "USDCHF": FXDeltaMethodTypes.SPOT_DELTA_PREM_ADJ,
    "EURCHF": FXDeltaMethodTypes.SPOT_DELTA_PREM_ADJ,
    "GBPUSD": FXDeltaMethodTypes.SPOT_DELTA,
    "EURGBP": FXDeltaMethodTypes.SPOT_DELTA_PREM_ADJ,
    "AUDUSD": FXDeltaMethodTypes.SPOT_DELTA,
    "AUDJPY": FXDeltaMethodTypes.SPOT_DELTA_PREM_ADJ,
    "USDCAD": FXDeltaMethodTypes.SPOT_DELTA_PREM_ADJ,
    "USDBRL": FXDeltaMethodTypes.SPOT_DELTA_PREM_ADJ,
    "USDMXN": FXDeltaMethodTypes.SPOT_DELTA_PREM_ADJ,
}

########################################################################################


class FinFXRate:

    def __init__(self, ccy1, ccy2, rate):

        if ccy1 in ccyNames:
            self.ccy1 = ccy1
        else:
            raise FinError("Unknown currency code:" + ccy1)

        if ccy2 in ccyNames:
            self.ccy2 = ccy2
        else:
            raise FinError("Unknown currency code:" + ccy2)

        self.ccy2 = ccy2
        self.rate = rate


########################################################################################
