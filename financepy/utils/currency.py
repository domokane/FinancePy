##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


########################################################################################

from .error import FinError
from .calendar import CalendarTypes

g_ccy_names = {
    "AED": ("UNITED ARAB EMIRATES", "DIRHAM"),
    "AUD": ("AUSTRALIA", "DOLLAR"),
    "BHD": ("BAHRAIN", "DINAR"),
    "BRL": ("BRAZIL", "REAL"),
    "CAD": ("CANADA", "DOLLAR"),
    "CHF": ("SWITZERLAND", "FRANC"),
    "CLP": ("CHILE", "PESO"),
    "CNH": ("OFFSHORE CHINESE RENMINBI", "YUAN"),
    "CNY": ("CHINA", "RENMINBI"),
    "COP": ("COLOMBIA", "PESO"),
    "CZK": ("CZECH REPUBLIC", "KORUNA"),
    "DKK": ("DENMARK", "KRONE"),
    "EUR": ("EUROZONE", "EURO"),
    "GBP": ("UK", "POUND"),
    "HKD": ("HONG KONG", "DOLLAR"),
    "HUF": ("HUNGARY", "FORINT"),
    "IDR": ("INDONESIA", "RUPIAH"),
    "INR": ("INDIA", "RUPEE"),
    "ILS": ("ISRAEL", "SHEKEL"),
    "JOD": ("JORDAN", "DINAR"),
    "JPY": ("JAPAN", "YEN"),
    "KRW": ("KOREAN", "WON"),
    "KWD": ("KUWAIT", "DINAR"),
    "MAD": ("MOROCCO", "DIRHAM"),
    "MYR": ("MALAYSIA", "RINGGIT"),
    "MXN": ("MEXICO", "PESO"),
    "NOK": ("NORWAY", "KRONE"),
    "NZD": ("NEW ZEALAND", "DOLLAR"),
    "OMR": ("OMAN", "RIAL"),
    "PHP": ("PHILIPPINES", "PESO"),
    "PLN": ("POLAND", "ZLOTY"),
    "QAR": ("QATAR", "RIYAL"),
    "RON": ("ROMANIA", "LEU"),
    "RUB": ("RUSSIA", "RUBLE"),
    "SAR": ("SAUDI ARABIA", "RIYAL"),
    "SEK": ("SWEDEN", "KRONA"),
    "SGD": ("SINGAPORE", "DOLLAR"),
    "THB": ("THAILAND", "BAHT"),
    "TND": ("TUNISIA", "DINAR"),
    "TRY": ("TURKEY", "LIRA"),
    "TWD": ("TAIWAN", "DOLLAR"),
    "USD": ("US", "DOLLAR"),
    "ZAR": ("SOUTH AFRICA", "RAND"),
}


CCY_CALENDARS = {
    "USD": CalendarTypes.UNITED_STATES,
    "EUR": CalendarTypes.TARGET,
    "GBP": CalendarTypes.UNITED_KINGDOM,
    "JPY": CalendarTypes.JAPAN,
    "CHF": CalendarTypes.SWITZERLAND,
    "CAD": CalendarTypes.CANADA,
    "AUD": CalendarTypes.AUSTRALIA,
    "NZD": CalendarTypes.NEW_ZEALAND,
}


def get_currency_calendar_type(currency: str) -> CalendarTypes:
    """Return the holiday calendar associated with a currency."""

    currency = currency.upper()

    try:
        return CCY_CALENDARS[currency]
    except KeyError:
        raise FinError(
            f"No holiday calendar defined for currency: {currency}"
        ) from None
